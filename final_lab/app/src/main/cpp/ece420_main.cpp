//
// Created by daran on 1/12/2017 to be used in ECE420 Sp17 for the first time.
// Modified by dwang49 on 1/1/2018 to adapt to Android 7.0 and Shield Tablet updates.
// Modified by rra2 and wonjoon2 on 04/15/2022 to adapt for final project
//

#include "ece420_main.h"
#include "ece420_lib.h"
#include "kiss_fft/kiss_fft.h"
#include <vector>

// JNI Function
extern "C" {
JNIEXPORT jintArray JNICALL
Java_com_ece420_lab4_MainActivity_getFreqUpdate(JNIEnv *env, jclass);
}

// Student Variables
#define ALPHA 0.97
#define FRAME_SIZE 1024 // 21.3333ms [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
#define F_S 48000
#define N_CEP 12
#define N_FFT 512
#define NUM_CLUSTERS 3
#define NUM_AUDIO_FRAMES 50
#define NUM_FILTER_BANK 20
#define F_BANK_SIZE N_FFT/2 + 1
#define NUM_CEP_COEF 12
#define OVERLAP_RATIO 0.5
#define NUM_OVERLAPPING_FRAMES int( NUM_AUDIO_FRAMES / OVERLAP_RATIO ) - 1 //500 * 2 = 1000 frames
#define NUM_OVERLAPPING_FRAME_NUM_CEP_COEF 1188
#define STRIDE_SIZE int( FRAME_SIZE * OVERLAP_RATIO ) //10ms --> half of the frame overlaps with the next frame
#define VOICED_THRESHOLD 100  // Find your own threshold

//Global Variables
int prediction_labels[NUM_OVERLAPPING_FRAMES] = { 0 };
float audioBuffer[FRAME_SIZE * NUM_AUDIO_FRAMES]; //20ms *
int audiobuffer_frameIdx = 0;
int audiobuffer_isfull = 0;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Adapted from: https://www.geeksforgeeks.org/discrete-cosine-transform-algorithm-program/
// Sourced date: 04/15/2022
// CPP program to perform discrete cosine transform

#define M NUM_OVERLAPPING_FRAMES
#define N NUM_FILTER_BANK

// Function to find discrete cosine transform and print it
//template<class T>
void dctTransform(float (&matrix)[M][N], float (&dct)[M][N])
{
    int i, j, k, l;

    float ci, cj, dct1, sum;
    const int m = M, n = N;

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {

            // ci and cj depends on frequency as well as
            // number of row and columns of specified matrix
            if (i == 0)
                ci = (float) 1.0 / sqrt(m);
            else
                ci = (float) sqrt(2.0) / sqrt(m);
            if (j == 0)
                cj = (float) 1.0 / sqrt(n);
            else
                cj = (float) sqrt(2.0) / sqrt(n);

            // sum will temporarily store the sum of
            // cosine signals
            sum = 0;
            for (k = 0; k < m; k++) {
                for (l = 0; l < n; l++) {
                    dct1 = matrix[k][l] *
                           cos((2 * k + 1) * i * M_PI / (2 * m)) *
                           cos((2 * l + 1) * j * M_PI / (2 * n));
                    sum += dct1;
                }
            }
            dct[i][j] = ci * cj * sum;
        }
    }
}
//This code is contributed by SoumikMondal

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Adapted from: https://www.geeksforgeeks.org/c-program-multiply-two-matrices/
// Sourced date: 04/15/2022
 /*
 * This C++ program can multiply any two square or rectangular matrices.
 * The below program multiplies two square matrices of size 4 * 4.
 * There is also an example of a rectangular matrix for the same code (commented below).
 * We can change the Matrix value with the number of rows and columns (from MACROs) for Matrix-1
 * and Matrix-2 for different dimensions.
 */

 /*
 * Note:  i- The number of columns in Matrix-1 must be equal to the number of rows in Matrix-2.
 *       ii- Output of multiplication of Matrix-1 and Matrix-2, results with equal to the number
 *           of rows of Matrix-1 and the number of columns of Matrix-2 i.e. rslt[R1][C2].
 */

// Edit MACROs here, according to your Matrix Dimensions for mat1[R1][C1] and mat2[R2][C2]
#define R1 NUM_OVERLAPPING_FRAMES // number of rows in Matrix-1
#define C1 F_BANK_SIZE            // number of columns in Matrix-1
#define R2 F_BANK_SIZE            // number of rows in Matrix-2
#define C2 NUM_FILTER_BANK        // number of columns in Matrix-2

template<class T>
void mulMat(T (&mat1)[R1][C1], T (&mat2)[R2][C2], T (&rslt)[R1][C2]) {

    for (int i = 0; i < R1; i++) {
        for (int j = 0; j < C2; j++) {
            rslt[i][j] = 0;
            for (int k = 0; k < R2; k++) {
                rslt[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
}

// This code is contributed by Manish Kumar (mkumar2789)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


// Helper Functions:
// Converts a frequency [Hz] into mel_frequency
float freq_to_mel(float frequency){
    return 2595 * log10(1 + (frequency/700));
}

float to_dB(float data){
    return 20 * log10(data);
}

void mel_to_freq(float * data, int mel_size){
//    float freq[mel_size];
    for (int i = 0; i < mel_size; i++) {
        float exponent = (data[i]/2595) - 1;
        data[i] = 700 * (pow(10, exponent));
    }
}

void noise_removal(float (&audioIn)[FRAME_SIZE * NUM_AUDIO_FRAMES], float threshold){
    // Voice Detection

    for(int i = 0; i < FRAME_SIZE * NUM_AUDIO_FRAMES; i++){
        if(audioIn[i] * audioIn[i] < threshold){
            audioIn[i] = 1;
        }
    }
}

void pre_emph(float (&data)[FRAME_SIZE * NUM_AUDIO_FRAMES], float alpha){
    for(int i = 1; i < FRAME_SIZE * NUM_AUDIO_FRAMES; i++){
        data[i] = data[i] - alpha * data[i - 1];
    }
}

void framing_windowing_power_spectrum(float (&audioIn)[FRAME_SIZE * NUM_AUDIO_FRAMES], float (&framed_audio)[NUM_OVERLAPPING_FRAMES][FRAME_SIZE], float (&framed_FFT)[NUM_OVERLAPPING_FRAMES][F_BANK_SIZE]){

    for(int frameIdx = 0; frameIdx < NUM_OVERLAPPING_FRAMES; frameIdx++){
        kiss_fft_cpx original_time[FRAME_SIZE];
        kiss_fft_cpx original_FREQUENCY[N_FFT]; // TODO: Potenital bug --> replace FRAME_SIZE w/ NFFT = 512 or 256
        kiss_fft_cfg FFT = kiss_fft_alloc(FRAME_SIZE, 0, NULL,NULL);

        for(int dataIdx = 0; dataIdx < FRAME_SIZE; dataIdx++){
            framed_audio[frameIdx][dataIdx] = audioIn[(frameIdx * STRIDE_SIZE) + dataIdx];
            framed_audio[frameIdx][dataIdx] *= 0.54 - ( 0.46 * cos((2 * M_PI * dataIdx)/(FRAME_SIZE - 1)));

            original_time[dataIdx].r = framed_audio[frameIdx][dataIdx];
            original_time[dataIdx].i = 0;
        }

        kiss_fft(FFT, original_time, original_FREQUENCY);
        for(int FFTdataIdx = 0; FFTdataIdx < F_BANK_SIZE; FFTdataIdx++){
            framed_FFT[frameIdx][FFTdataIdx] = original_FREQUENCY[FFTdataIdx].r * original_FREQUENCY[FFTdataIdx].r;
            framed_FFT[frameIdx][FFTdataIdx] += original_FREQUENCY[FFTdataIdx].i * original_FREQUENCY[FFTdataIdx].i;
            framed_FFT[frameIdx][FFTdataIdx] /= N_FFT;
        }
    }
}

void filterbank(float (&fbank)[F_BANK_SIZE][NUM_FILTER_BANK]){
    float lower_bound = 0;
    float upper_bound = freq_to_mel(float(F_S));
    float step = (upper_bound-lower_bound)/(NUM_FILTER_BANK + 2);
    int mel_size = NUM_FILTER_BANK + 2;

    float mel_scale[mel_size];
    for(int i = 0; i < mel_size; i++){
        mel_scale[i] = (step * i) + lower_bound;
    }

    mel_to_freq(mel_scale, mel_size);

    int bin_f[NUM_FILTER_BANK + 2] = {0};
    for(int i = 0; i < mel_size; i++){
        bin_f[i] = int(mel_scale[i] * int(N_FFT/2 + 1) / F_S);
    }

    for(int m = 1; m <= NUM_FILTER_BANK; m++){
        for(int k = 0; k < int(F_BANK_SIZE); k++){
            if(k < bin_f[m-1]){
                fbank[k][m-1] = 0;
            }
            if(k >= bin_f[m-1] && k < bin_f[m]){
                fbank[k][m-1] = (float)((k - bin_f[m-1])/(bin_f[m] - bin_f[m-1]));
            }
            if(k >= bin_f[m] && k < bin_f[m+1]){
                fbank[k][m-1] = (float)(bin_f[m+1] - k)/(bin_f[m+1] - bin_f[m]);
            }
            if(k >= bin_f[m+1]){
                fbank[k][m-1] = 0;
            }
        }
    }


}

void mfcc(float (&audioIn)[FRAME_SIZE * NUM_AUDIO_FRAMES], int audio_size, float (&coef)[NUM_OVERLAPPING_FRAMES][NUM_FILTER_BANK]){
    noise_removal(audioIn, VOICED_THRESHOLD);
    pre_emph(audioIn, ALPHA);

    float framed_audio[NUM_OVERLAPPING_FRAMES][FRAME_SIZE] = {{0}};
    float power_spectrum[NUM_OVERLAPPING_FRAMES][F_BANK_SIZE] = {{0}};
    framing_windowing_power_spectrum(audioIn, framed_audio, power_spectrum);

    float fbank[F_BANK_SIZE][NUM_FILTER_BANK] = {{0}};
    filterbank(fbank); //TODO: Verify if the bandpass filtering is correct

    //    pre_coef = power_spectrum @ filter_bank
    float pre_coef[NUM_OVERLAPPING_FRAMES][NUM_FILTER_BANK] = {{0}};
    mulMat(power_spectrum, fbank, pre_coef);

    //    Convert to decibel scale (dB)
    for(int i = 0; i < NUM_OVERLAPPING_FRAMES; i++){
        for(int j = 0; j < NUM_FILTER_BANK; j++){
            pre_coef[i][j] = to_dB(pre_coef[i][j]);
            if( isinf(pre_coef[i][j]) ){
                pre_coef[i][j] = -99999999.9;
            }
        }
    }

    //    DCT of audio frames
    dctTransform(pre_coef, coef);

    //    Normalize Audio frames
    for(int frameIdx = 0; frameIdx < NUM_OVERLAPPING_FRAMES; frameIdx++){
        float mean = 0;
        for(int coefIdx = 0; coefIdx < NUM_FILTER_BANK; coefIdx++){
            mean += coef[frameIdx][coefIdx];
        }
        mean /= NUM_FILTER_BANK;
        mean += 0.00000001; // 1e-8
        for(int coefIdx = 0; coefIdx < NUM_FILTER_BANK; coefIdx++){
            coef[frameIdx][coefIdx] -= mean;
        }
    }
}

class KMeans {
    private:
        int iterations;
        float centers[NUM_CLUSTERS][NUM_CEP_COEF];
        int labels[NUM_OVERLAPPING_FRAMES];
    public:

        // Computes Frobenius norm of the matrix
//        void norm(float (&data)[NUM_OVERLAPPING_FRAMES][NUM_CEP_COEF], float (&result)[NUM_OVERLAPPING_FRAMES]){
//            for(int frameIdx = 0; frameIdx < NUM_OVERLAPPING_FRAMES; frameIdx++){
//                for(int cepIdx = 0; cepIdx < NUM_CEP_COEF; cepIdx++){
//                    result[frameIdx] += abs(data[frameIdx][cepIdx]) * abs(data[frameIdx][cepIdx]);
//                }
//            }
//        }

        KMeans(){
            iterations = 100;
        }

        KMeans(int iter){
            iterations = iter;
        }

        void dist_calc(float (&vectors)[NUM_OVERLAPPING_FRAMES][NUM_CEP_COEF], float (&distance)[NUM_OVERLAPPING_FRAMES][NUM_CLUSTERS]){

            float dif[NUM_OVERLAPPING_FRAMES][NUM_CEP_COEF] = {{0}};

            for(int clustIdx = 0; clustIdx < NUM_CLUSTERS; clustIdx++){
                for(int frameIdx = 0; frameIdx < NUM_OVERLAPPING_FRAMES; frameIdx++){
                    for(int cepIdx = 0; cepIdx < NUM_CEP_COEF; cepIdx++){
                        dif[frameIdx][cepIdx] = vectors[frameIdx][cepIdx] - this->centers[clustIdx][cepIdx];
                    }
                    for(int cepIdx = 0; cepIdx < NUM_CEP_COEF; cepIdx++){
                        distance[frameIdx][clustIdx] += abs(dif[frameIdx][cepIdx]) * abs(dif[frameIdx][cepIdx]);
                    }
                }
            }

        }

        void center_init(float (&vectors)[NUM_OVERLAPPING_FRAMES][NUM_CEP_COEF]){
            //Randomize initial center audio vectors
            std::vector<int> frameIdxs(NUM_OVERLAPPING_FRAMES, 0);
            for(int idx = 0; idx < NUM_OVERLAPPING_FRAMES; idx++){
                frameIdxs[idx] = idx;
            }
            std::random_shuffle(frameIdxs.begin(), frameIdxs.end());

            int frameIdx = 0;
            for(int c = 0; c < NUM_CLUSTERS; c++){
                frameIdx = frameIdxs[c];
                for(int cepIdx = 0; cepIdx < NUM_CEP_COEF; cepIdx++) {
                    this->centers[c][cepIdx] = vectors[frameIdx][cepIdx];
                }
            }
        }

        void center_init_pretrained(){
            // Gender Classification
            float pretrained_centers[NUM_CLUSTERS][NUM_CEP_COEF] = {{-14.66965037,10.96730311, 0.61583515,-8.35761944, 2.57579316, -4.5830111, 12.16768616,-7.05242872, 0.80650813,-1.0964556,0.48433202, 3.39654526},
                    { -0.87859462,11.44264234,-19.41220519,22.64120235,6.07230454,-12.33452725,-5.71475363,-1.97208594,12.79446273,-2.00796186, -6.88701148, 1.04113228},
                    { 18.67582217 -20.49063213,10.98780174,-3.36683004,-6.86119169, 13.13409599, -11.59057677, 9.91637172,-8.74154015, 2.57148776, 3.56941206,-4.83114923}};
            for(int c = 0; c < NUM_CLUSTERS; c++) {
                for (int cepIdx = 0; cepIdx < NUM_CEP_COEF; cepIdx++) {
                    this->centers[c][cepIdx] = pretrained_centers[c][cepIdx];
                }
            }
        }

        void update_centers(float (&vectors)[NUM_OVERLAPPING_FRAMES][NUM_CEP_COEF]){
            for(int c = 0; c < NUM_CLUSTERS; c++) {
                std::vector<std::vector<float>> cluster_vectors;
                for(int frameIdx = 0; frameIdx < NUM_OVERLAPPING_FRAMES; frameIdx++) {
                    if(this->labels[frameIdx] == c) {
                        std::vector<float> frame;
                        for(int i = 0; i < NUM_CEP_COEF; i++){
                            frame.push_back(vectors[frameIdx][i]);
                        }
                        cluster_vectors.push_back(frame);
                    }
                }
                float mean[NUM_CEP_COEF] = { 0 };
                for(std::vector<float>frame : cluster_vectors){
                    for(int coef = 0; coef < NUM_CEP_COEF; coef++){
                        mean[coef] += frame[coef];
                    }
                }
                for(int coef = 0; coef < NUM_CEP_COEF; coef++){
                    mean[coef] /= NUM_OVERLAPPING_FRAMES;
                    this->centers[c][coef] = mean[coef];
                }
            }
        }
        void find_closest(float (&distance)[NUM_OVERLAPPING_FRAMES][NUM_CLUSTERS], int (&result)[NUM_OVERLAPPING_FRAMES]){

            for(int frameIdx = 0; frameIdx < NUM_OVERLAPPING_FRAMES; frameIdx++) {
                float min_distance = FLT_MAX;
                int min_distanceIdx = 0;
                for (int c = 0; c < NUM_CLUSTERS; c++) {
                    if(distance[frameIdx][c] < min_distance){
                        min_distance = distance[frameIdx][c];
                        min_distanceIdx = c;
                    }
                }
                result[frameIdx] = min_distanceIdx;
            }
        }
        void fit_data(float (&vectors)[NUM_OVERLAPPING_FRAMES][NUM_CEP_COEF]){
            float dist[NUM_OVERLAPPING_FRAMES][NUM_CLUSTERS];

            center_init(vectors);

            for(int iter = 0; iter < this->iterations; iter++){
                dist_calc(vectors, dist);
                find_closest(dist, this->labels);
                update_centers(vectors);
            }
        }
        void predict(float (&vectors)[NUM_OVERLAPPING_FRAMES][NUM_CEP_COEF], int (&predictions)[NUM_OVERLAPPING_FRAMES]){
            float dist[NUM_OVERLAPPING_FRAMES][NUM_CLUSTERS];

            dist_calc(vectors, dist);
            find_closest(dist, predictions);
        }
};

void speaker_identification(){
    float mel_coef[NUM_OVERLAPPING_FRAMES][NUM_FILTER_BANK] = {{0}};
    float mel_coef_reduced[NUM_OVERLAPPING_FRAMES][NUM_CEP_COEF] = {{0}};
    float mel_coef_dump[NUM_OVERLAPPING_FRAME_NUM_CEP_COEF] = {0};

    mfcc(audioBuffer, NUM_AUDIO_FRAMES * FRAME_SIZE, mel_coef);

    for(int f = 0; f < NUM_OVERLAPPING_FRAMES; f++) {
        for (int c = 0; c < NUM_CEP_COEF; c++) {
            if(f <= 1){
                if(mel_coef[2][c + 1] >= 50.0){
                    mel_coef[2][c + 1] = 50.0;
                }
                if(mel_coef[2][c + 1] <= -50.0){
                    mel_coef[2][c + 1] = -50.0;
                }

                mel_coef_reduced[f][c] = mel_coef[2][c + 1];
                mel_coef_dump[(f * NUM_CEP_COEF) + c] = mel_coef[2][c + 1];
                continue;
            }
            if(mel_coef[f][c + 1] >= 50.0){
                mel_coef[f][c + 1] = 50.0;
            }
            if(mel_coef[f][c + 1] <= -50.0){
                mel_coef[f][c + 1] = -50.0;
            }
            mel_coef_reduced[f][c] = mel_coef[f][c + 1];
            mel_coef_dump[(f * NUM_CEP_COEF) + c] = mel_coef[f][c + 1];
        }
    }
//    int prediction_labels[NUM_OVERLAPPING_FRAMES] = {0};
    KMeans kmeans = KMeans(100);

//    kmeans.fit_data(mel_coef_reduced);
    kmeans.center_init_pretrained();
    kmeans.predict(mel_coef_reduced, prediction_labels);
}

void ece420ProcessFrame(sample_buf *dataBuf) {
    // Keep in mind, we only have 20ms to process each buffer!
    struct timeval start;
    struct timeval end;
    gettimeofday(&start, NULL);

    // Data is encoded in signed PCM-16, little-endian, mono
    float bufferIn[FRAME_SIZE];
    for (int i = 0; i < FRAME_SIZE; i++) {
        int16_t val = ((uint16_t) dataBuf->buf_[2 * i]) | (((uint16_t) dataBuf->buf_[2 * i + 1]) << 8);
        bufferIn[i] = (float) val;
        if(audiobuffer_isfull == 0){
            audioBuffer[(audiobuffer_frameIdx*FRAME_SIZE) + i] = (float) val;
        }
    }

    audiobuffer_frameIdx++;
    if(audiobuffer_frameIdx >= NUM_AUDIO_FRAMES){
        audiobuffer_isfull = 1;
    }

    gettimeofday(&end, NULL);
    LOGD("Time delay: %ld us",  ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)));
}

extern "C" JNIEXPORT jintArray JNICALL
Java_com_ece420_lab4_MainActivity_getFreqUpdate(JNIEnv *env, jclass) {
//  Modified from: https://stackoverflow.com/questions/1610045/how-to-return-an-array-from-jni-to-java
//  Date Accessed:
    jintArray result;
    result = (*env).NewIntArray(NUM_OVERLAPPING_FRAMES);
    jint fill[NUM_OVERLAPPING_FRAMES];
    for (int i = 0; i < NUM_OVERLAPPING_FRAMES; i++) {
        fill[i] = prediction_labels[i];
    }

    (*env).SetIntArrayRegion(result, 0, NUM_OVERLAPPING_FRAMES, fill);
    return result;
}