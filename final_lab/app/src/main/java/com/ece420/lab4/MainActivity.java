/*
 * Copyright 2015 The Android Open Source Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.ece420.lab4;

import android.app.Activity;
import android.content.Context;
import android.content.pm.ActivityInfo;
import android.content.pm.PackageManager;
import android.Manifest;
import android.media.AudioFormat;
import android.media.AudioManager;
import android.media.AudioRecord;
import android.os.AsyncTask;
import android.os.Bundle;
import android.os.Handler;
import android.support.annotation.NonNull;
import android.support.v4.app.ActivityCompat;
import android.view.Menu;
import android.view.MenuItem;
import android.view.View;
import android.view.WindowManager;
import android.widget.Button;
import android.widget.TextView;
import android.widget.Toast;
import android.widget.ProgressBar;

import java.io.File;

import java.util.Arrays;
import java.util.Timer;
import java.util.TimerTask;

//import androidx.appcompat.app.AppCompatActivity;

import com.jjoe64.graphview.GraphView;
import com.jjoe64.graphview.series.DataPoint;
import com.jjoe64.graphview.series.LineGraphSeries;


public class MainActivity extends Activity
        implements ActivityCompat.OnRequestPermissionsResultCallback {

    // UI Variables
    Button   controlButton;
    TextView statusView;
    ProgressBar simpleProgressBar;
    int progressValue;
    static TextView freq_view;
    String  nativeSampleRate;
    String  nativeSampleBufSize;
    boolean supportRecording;
    GraphView graphView;
    // Static Values
    Boolean isPlaying = false;
    private static final int AUDIO_ECHO_REQUEST = 0;
    private static final int FRAME_SIZE = 1024;
    private static final int PROGRESS_MAX = 100;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        getWindow().addFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON);
        setContentView(R.layout.activity_main);
        super.setRequestedOrientation (ActivityInfo.SCREEN_ORIENTATION_PORTRAIT);

        // Google NDK Stuff
        controlButton = (Button)findViewById((R.id.capture_control_button));
        statusView = (TextView)findViewById(R.id.statusView);
        simpleProgressBar=(ProgressBar)findViewById(R.id.simpleProgressBar); // initiate the progress bar
        simpleProgressBar.setMax(PROGRESS_MAX);
        progressValue = 0;
        simpleProgressBar.setProgress(progressValue);
        queryNativeAudioParameters();
        // initialize native audio system
        updateNativeAudioUI();
        if (supportRecording) {
            createSLEngine(Integer.parseInt(nativeSampleRate), FRAME_SIZE);
        }
        // on below line we are initializing our graph view.
//        graphView = (GraphView) findViewById(R.id.idGraphView);
//
//        // on below line we are adding data to our graph view.
//        LineGraphSeries<DataPoint> series = new LineGraphSeries<DataPoint>(new DataPoint[]{
//                // on below line we are adding
//                // each point on our x and y axis.
////                new DataPoint(0, 1),
////                new DataPoint(1, 3),
////                new DataPoint(2, 4),
////                new DataPoint(3, 9),
////                new DataPoint(4, 6),
////                new DataPoint(5, 3),
////                new DataPoint(6, 6),
////                new DataPoint(7, 1),
////                new DataPoint(8, 2)
//        });
//
//        // after adding data to our line graph series.
//        // on below line we are setting
//        // title for our graph view.
//        graphView.setTitle("My Graph View");
//
//        // on below line we are setting
//        // text color to our graph view.
////        graphView.setTitleColor(R.color.purple_200);
//
//        // on below line we are setting
//        // our title text size.
//        graphView.setTitleTextSize(18);
//
//        // on below line we are adding
//        // data series to our graph view.
//        graphView.addSeries(series);

        // Setup UI
        freq_view = (TextView)findViewById(R.id.textFrequency);
        initializeFreqTextBackgroundTask(100);
    }
    @Override
    protected void onDestroy() {
        if (supportRecording) {
            if (isPlaying) {
                stopPlay();
            }
            deleteSLEngine();
            isPlaying = false;
        }
        super.onDestroy();
    }

    @Override
    public boolean onCreateOptionsMenu(Menu menu) {
        // Inflate the menu; this adds items to the action bar if it is present.
        getMenuInflater().inflate(R.menu.menu_main, menu);
        return true;
    }

    @Override
    public boolean onOptionsItemSelected(MenuItem item) {
        // Handle action bar item clicks here. The action bar will
        // automatically handle clicks on the Home/Up button, so long
        // as you specify a parent activity in AndroidManifest.xml.
        int id = item.getItemId();

        //noinspection SimplifiableIfStatement
        if (id == R.id.action_settings) {
            return true;
        }

        return super.onOptionsItemSelected(item);
    }

    private void startEcho() {
        if(!supportRecording){
            return;
        }
        if (!isPlaying) {
            if(!createSLBufferQueueAudioPlayer()) {
                statusView.setText(getString(R.string.error_player));
                return;
            }
            if(!createAudioRecorder()) {
                deleteSLBufferQueueAudioPlayer();
                statusView.setText(getString(R.string.error_recorder));
                return;
            }
            startPlay();   // this must include startRecording()
            statusView.setText(getString(R.string.status_echoing));
        } else {
            stopPlay();  //this must include stopRecording()
            progressValue = 0;
            simpleProgressBar.setProgress(progressValue);
            updateNativeAudioUI();
            deleteAudioRecorder();
            deleteSLBufferQueueAudioPlayer();
        }
        isPlaying = !isPlaying;
        controlButton.setText(getString((isPlaying == true) ?
                R.string.StopEcho: R.string.StartEcho));
    }

    public void onEchoClick(View view) {
        if (ActivityCompat.checkSelfPermission(this, Manifest.permission.RECORD_AUDIO) !=
                PackageManager.PERMISSION_GRANTED) {
            statusView.setText(getString(R.string.status_record_perm));
            ActivityCompat.requestPermissions(
                    this,
                    new String[] { Manifest.permission.RECORD_AUDIO },
                    AUDIO_ECHO_REQUEST);
            return;
        }
        startEcho();
    }

    public void getLowLatencyParameters(View view) {
        updateNativeAudioUI();
        return;
    }

    private void queryNativeAudioParameters() {
        AudioManager myAudioMgr = (AudioManager) getSystemService(Context.AUDIO_SERVICE);
        nativeSampleRate  =  myAudioMgr.getProperty(AudioManager.PROPERTY_OUTPUT_SAMPLE_RATE);
        nativeSampleBufSize =myAudioMgr.getProperty(AudioManager.PROPERTY_OUTPUT_FRAMES_PER_BUFFER);
        int recBufSize = AudioRecord.getMinBufferSize(
                Integer.parseInt(nativeSampleRate),
                AudioFormat.CHANNEL_IN_MONO,
                AudioFormat.ENCODING_PCM_16BIT);
        supportRecording = true;
        if (recBufSize == AudioRecord.ERROR ||
                recBufSize == AudioRecord.ERROR_BAD_VALUE) {
            supportRecording = false;
        }
    }
    private void updateNativeAudioUI() {
        if (!supportRecording) {
            statusView.setText(getString(R.string.error_no_mic));
            controlButton.setEnabled(false);
            return;
        }

        statusView.setText("nativeSampleRate    = " + nativeSampleRate + "\n" +
                "nativeSampleBufSize = " + nativeSampleBufSize + "\n");

    }
    @Override
    public void onRequestPermissionsResult(int requestCode, @NonNull String[] permissions,
                                           @NonNull int[] grantResults) {
        /*
         * if any permission failed, the sample could not play
         */
        if (AUDIO_ECHO_REQUEST != requestCode) {
            super.onRequestPermissionsResult(requestCode, permissions, grantResults);
            return;
        }

        if (grantResults.length != 1  ||
                grantResults[0] != PackageManager.PERMISSION_GRANTED) {
            /*
             * When user denied permission, throw a Toast to prompt that RECORD_AUDIO
             * is necessary; also display the status on UI
             * Then application goes back to the original state: it behaves as if the button
             * was not clicked. The assumption is that user will re-click the "start" button
             * (to retry), or shutdown the app in normal way.
             */
            statusView.setText(getString(R.string.error_no_permission));
            Toast.makeText(getApplicationContext(),
                    getString(R.string.prompt_permission),
                    Toast.LENGTH_SHORT).show();
            return;
        }

        /*
         * When permissions are granted, we prompt the user the status. User would
         * re-try the "start" button to perform the normal operation. This saves us the extra
         * logic in code for async processing of the button listener.
         */
        statusView.setText("RECORD_AUDIO permission granted, touch " +
                getString(R.string.StartEcho) + " to begin");

        // The callback runs on app's thread, so we are safe to resume the action
        startEcho();
    }

    // All this does is calls the UpdateStftTask at a fixed interval
    // http://stackoverflow.com/questions/6531950/how-to-execute-async-task-repeatedly-after-fixed-time-intervals
    public void initializeFreqTextBackgroundTask(int timeInMs) {
        final Handler handler = new Handler();
        Timer timer = new Timer();
        TimerTask doAsynchronousTask = new TimerTask() {
            @Override
            public void run() {
                handler.post(new Runnable() {
                    public void run() {
                        try {
                            UpdateFreqTextTask performFreqTextUpdate = new UpdateFreqTextTask();
                            performFreqTextUpdate.execute();
                        } catch (Exception e) {
                            // TODO Auto-generated catch block
                        }
                    }
                });
            }
        };
        timer.schedule(doAsynchronousTask, 0, timeInMs); // execute every 100 ms
    }

    // UI update
    private class UpdateFreqTextTask extends AsyncTask<Void, int[], Void> {
        @Override
        protected Void doInBackground(Void... params) {

            // Update screen, needs to be done on UI thread
            publishProgress(getFreqUpdate());

            return null;
        }

        protected void onProgressUpdate(int[]... newFreq) {
            if(isPlaying) {
                progressValue++;
                if (progressValue >= PROGRESS_MAX) {
                    progressValue = PROGRESS_MAX;
                }
                simpleProgressBar.setProgress(progressValue);
            }

//            LineGraphSeries<DataPoint> series = new LineGraphSeries<DataPoint>(new DataPoint[newFreq[0].length]);

            // after adding data to our line graph series.
            // on below line we are setting
            // title for our graph view.
//            graphView.setTitle("My Graph View");

            // on below line we are setting
            // text color to our graph view.
//        graphView.setTitleColor(R.color.purple_200);

            String text, output = "";
            int data;
            for (int i = 0; i < newFreq[0].length; i++){
                data = newFreq[0][i];
//                DataPoint datap = new DataPoint(i, data);
//                series.;
                switch (data) {
                    case 0:
                        text = "Background";
                        break;
                    case 1:
                        text = "Male";
                        break;
                    case 2:
                        text = "Female";
                        break;
                    default:
                        text = "";
                        break;
                }
                output += text + ", ";
            }

            freq_view.setText(output);
            // on below line we are setting
//            // our title text size.
//            graphView.setTitleTextSize(18);
//
//            // on below line we are adding
//            // data series to our graph view.
//            graphView.addSeries(series);

//            if (newFreq[0] > 0) {
//                freq_view.setText(Long.toString(newFreq.intValue()) + " Hz");
//            } else {
//                freq_view.setText("Unvoiced");
//            }
        }
    }

    /*
     * Loading our Libs
     */
    static {
        System.loadLibrary("echo");
    }

    /*
     * jni function implementations...
     */
    public static native void createSLEngine(int rate, int framesPerBuf);
    public static native void deleteSLEngine();

    public static native boolean createSLBufferQueueAudioPlayer();
    public static native void deleteSLBufferQueueAudioPlayer();

    public static native boolean createAudioRecorder();
    public static native void deleteAudioRecorder();
    public static native void startPlay();
    public static native void stopPlay();

    public static native int[] getFreqUpdate();
}
