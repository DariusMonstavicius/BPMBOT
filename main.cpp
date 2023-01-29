/** @file paex_record.c
	@ingroup examples_src
	@brief Record input into an array; Save array to a file; Playback recorded data.
	@author Phil Burk  http://www.softsynth.com
*/
/*
 * $Id$
 *
 * This program uses the PortAudio Portable Audio Library.
 * For more information see: http://www.portaudio.com
 * Copyright (c) 1999-2000 Ross Bencina and Phil Burk
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files
 * (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge,
 * publish, distribute, sublicense, and/or sell copies of the Software,
 * and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
 * ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

 /*
  * The text above constitutes the entire PortAudio license; however,
  * the PortAudio community also makes the following non-binding requests:
  *
  * Any person wishing to distribute modifications to the Software is
  * requested to send the modifications to the original developer so that
  * they can be incorporated into the canonical version. It is also
  * requested that these non-binding requests be included along with the
  * license above.
  */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <algorithm>
#include <Windows.h>

using namespace std;
#include "portaudio.h"
#include "plot.h"



  /* #define SAMPLE_RATE  (17932) // Test failure to open with this value. */
#define SAMPLE_RATE  (44100)
#define FRAMES_PER_BUFFER (512)
#define NUM_SECONDS     (5)
#define NUM_CHANNELS    (1)
/* #define DITHER_FLAG     (paDitherOff) */
#define DITHER_FLAG     (0) /**/
/** Set to 1 if you want to capture the recording to a file. */
#define WRITE_TO_FILE   (0)

/* Select sample format. */
#if 1
#define PA_SAMPLE_TYPE  paFloat32
typedef float SAMPLE;
#define SAMPLE_SILENCE  (0.0f)
#define PRINTF_S_FORMAT "%.8f"
#elif 1
#define PA_SAMPLE_TYPE  paInt16
typedef short SAMPLE;
#define SAMPLE_SILENCE  (0)
#define PRINTF_S_FORMAT "%d"
#elif 0
#define PA_SAMPLE_TYPE  paInt8
typedef char SAMPLE;
#define SAMPLE_SILENCE  (0)
#define PRINTF_S_FORMAT "%d"
#else
#define PA_SAMPLE_TYPE  paUInt8
typedef unsigned char SAMPLE;
#define SAMPLE_SILENCE  (128)
#define PRINTF_S_FORMAT "%d"
#endif

typedef struct
{
	int          numbytes;
	int          isrecording; //callback will turn this on at start of recording and off at end of recording
	int          isprocessing;  //callback will wait for this to clear before recording new data
	int          frameIndex;  /* Index into sample array. */
	int          maxFrameIndex;
	SAMPLE      *recordedSamples;
}
paTestData;

/* This routine will be called by the PortAudio engine when audio is needed.
** It may be called at interrupt level on some machines so don't do anything
** that could mess up the system like calling malloc() or free().
*/
static int recordCallback(const void *inputBuffer, void *outputBuffer,
	unsigned long framesPerBuffer,
	const PaStreamCallbackTimeInfo* timeInfo,
	PaStreamCallbackFlags statusFlags,
	void *userData)
{
	paTestData *data = (paTestData*)userData;
	const SAMPLE *rptr = (const SAMPLE*)inputBuffer;
	SAMPLE *wptr = &data->recordedSamples[data->frameIndex * NUM_CHANNELS];
	
	long framesToCalc;
	long i;
	int finished;
	unsigned long framesLeft = data->maxFrameIndex - data->frameIndex;

	(void)outputBuffer; /* Prevent unused variable warnings. */
	(void)timeInfo;
	(void)statusFlags;
	(void)userData;
	while (data->isprocessing) {
		;  //wait for data to stop processing
	}
	data->isrecording = 1;
	if (framesLeft < framesPerBuffer)
	{
		framesToCalc = framesLeft;
		finished = paComplete;
	}
	else
	{
		framesToCalc = framesPerBuffer;
		finished = paContinue;
	}

	if (inputBuffer == NULL)
	{
		for (i = 0; i < framesToCalc; i++)
		{
			*wptr++ = SAMPLE_SILENCE;  /* left */
			if (NUM_CHANNELS == 2) *wptr++ = SAMPLE_SILENCE;  /* right */
		}
	}
	else
	{
		for (i = 0; i < framesToCalc; i++)
		{
			*wptr++ = *rptr++;  /* left */
			if (NUM_CHANNELS == 2) *wptr++ = *rptr++;  /* right */
		}
	}
	data->frameIndex += framesToCalc;
	data->isrecording = 0;

	//if we are finished with original buffer, set frame back to zero to start recording again at beginning.
	if (finished == paComplete) {
		data->frameIndex = 0;
	}

	return paContinue; // finished;
}

/* This routine will be called by the PortAudio engine when audio is needed.
** It may be called at interrupt level on some machines so don't do anything
** that could mess up the system like calling malloc() or free().
*/
static int playCallback(const void *inputBuffer, void *outputBuffer,
	unsigned long framesPerBuffer,
	const PaStreamCallbackTimeInfo* timeInfo,
	PaStreamCallbackFlags statusFlags,
	void *userData)
{
	paTestData *data = (paTestData*)userData;
	SAMPLE *rptr = &data->recordedSamples[data->frameIndex * NUM_CHANNELS];
	SAMPLE *wptr = (SAMPLE*)outputBuffer;
	unsigned int i;
	int finished;
	unsigned int framesLeft = data->maxFrameIndex - data->frameIndex;

	(void)inputBuffer; /* Prevent unused variable warnings. */
	(void)timeInfo;
	(void)statusFlags;
	(void)userData;

	if (framesLeft < framesPerBuffer)
	{
		/* final buffer... */
		for (i = 0; i < framesLeft; i++)
		{
			*wptr++ = *rptr++;  /* left */
			if (NUM_CHANNELS == 2) *wptr++ = *rptr++;  /* right */
		}
		for (; i < framesPerBuffer; i++)
		{
			*wptr++ = 0;  /* left */
			if (NUM_CHANNELS == 2) *wptr++ = 0;  /* right */
		}
		data->frameIndex += framesLeft;
		finished = paComplete;
	}
	else
	{
		for (i = 0; i < framesPerBuffer; i++)
		{
			*wptr++ = *rptr++;  /* left */
			if (NUM_CHANNELS == 2) *wptr++ = *rptr++;  /* right */
		}
		data->frameIndex += framesPerBuffer;
		finished = paContinue;
	}
	return finished;
}

//compute the absolute value of a number
SAMPLE  sampleabs(SAMPLE sample) {
	if (sample < 0) {
		sample = ((SAMPLE)(-1.0)) * sample;
	}
	return sample;
}

//compute the absolute value and average of an array of data
SAMPLE ComputeAbsAverageOfArray(SAMPLE *data, int data_len) {
	SAMPLE ave = (SAMPLE) 0.0;

	int i;
	for (i = 0; i < data_len; i++) {
		ave = ave + sampleabs(data[i]);
	}
	ave = ave / (SAMPLE)data_len;
	return ave;
}

//compute the absolute value and maximum of an array of data
SAMPLE ComputeAbsMaxOfArray(SAMPLE *data, int data_len, int *max_index) {
	SAMPLE max = (SAMPLE) 0.0;
	*max_index = 0;

	int i;
	for (i = 0; i < data_len; i++) {
		if (sampleabs(data[i]) > max) {
			max = sampleabs(data[i]);
			*max_index = i;
		}
	}

	return max;
}

//Computes the moving average of an array of data
void AverageFilter(SAMPLE *data, int data_len, int filterwidth) {
	int i;
	SAMPLE *tempdata = (SAMPLE *)malloc(data_len * sizeof(SAMPLE));
	memcpy(tempdata, data, data_len * sizeof(SAMPLE));
	for (i = 0; i < filterwidth; i++) {
		tempdata[i] = ComputeAbsAverageOfArray(data + i, filterwidth);
	}
	for (i = filterwidth; i < data_len - filterwidth; i++) {
		tempdata[i] = ComputeAbsAverageOfArray(data + i - filterwidth, filterwidth*2);
	}
	for (i = data_len - filterwidth; i < data_len; i++) {
		tempdata[i] = ComputeAbsAverageOfArray(data + i - filterwidth, filterwidth);
	}
	memcpy(data, tempdata, data_len * sizeof(SAMPLE));
	free(tempdata);
}

//Find the Indices (x values) of all samples above a value
std::vector<int> FindIndicesAboveValue(SAMPLE *data, int data_len, SAMPLE value) {
	std::vector<int> indices;

	int i;
	for (i = 0; i < data_len; i++) {
		if (sampleabs(data[i]) >= value) {
			indices.push_back(i);
		}
	}

	return indices;
}

//This function didn't work
std::vector<int> FindMaxIndices_old(std::vector<int> neighbors, SAMPLE *data, int data_len) {
	std::vector<int> indices;
	int i, j;
	SAMPLE max = (SAMPLE)0.0;

	for (i = 0, j = 0; i < neighbors.size(); i++) {
		j = i;
		while ((j < neighbors.size() - 1) && (neighbors[j + 1] == neighbors[j] + 1)) {
			j++;  //advance j to the next clump of neighbors
		}
		if (j != i) {
			indices.push_back(neighbors[i + (j - i) / 2]);
			i = j;
		}
	}

	return indices;
}

//Find the indices (x values) of the maximum
std::vector<int> FindMaxIndices(std::vector<int> neighbors, SAMPLE* data, int data_len) {
	std::vector<int> indices;
	int i, j;
	SAMPLE max = (SAMPLE)0.0;

	for (i = 0, j = 0; i < neighbors.size(); i++) {
		j = i;
		while ((j < neighbors.size() - 1) && (neighbors[j + 1] == neighbors[j] + 1)) {
			j++;  //advance j to the next clump of neighbors
		}
		if (j != i) {
			indices.push_back(neighbors[i + (j - i) / 2]);
			//int maxindex = i;
			//float max = ComputeAbsMaxOfArray(data + i, j - i, &maxindex);
			//indices.push_back(neighbors[maxindex]);
			i = j;
		}
	}

	return indices;
}

//find the distances between consecutive entries in a vector
std::vector<int> FindDistances(std::vector<int> maxindices) {
	std::vector<int> distances;
	int i;
	for (i = 0; i < maxindices.size() - 1; i++) {
		distances.push_back(maxindices[i + 1] - maxindices[i]);
	}
	return distances;
}

//compute the average of an array
SAMPLE FindAverageDistance(std::vector<int> distances) {
	SAMPLE ave = (SAMPLE) 0.0;
	int i;
	for (i = 0; i < distances.size(); i++) {
		ave += distances[i];
	}
	ave = ave / distances.size();
	return ave;
}

//coimpute the median of an array
int FindMedianDistance(std::vector<int> distances) {
	int median = 0;
	size_t size = distances.size();
	if (size > 0){
	  std::sort(distances.begin(), distances.end());
	  if (size % 2) {
		  median = distances[size / 2];  //odd case
	  } else {
		  median = (distances[size / 2-1] + distances[size / 2]) / 2;  //even case
	  }
	}
	return median;
}

//Function that cleans up the distances by dropping distances that are far away from the average
std::vector<int> CleanUpDistances(std::vector<int> distances, SAMPLE ave) {
	//only use distances that are near the average over the sample
	std::vector<int> cleandistances;
	int i;
	for (i = 0; i < distances.size(); i++) {
		if (sampleabs(distances[i] - ave) < ave / (SAMPLE) 4.0) {
			cleandistances.push_back(distances[i]);
		}
	}
	return cleandistances;
}

//main function that computes the BPM
float ComputeBPM(paTestData *data) {
	float distancebetweenpeaks, bpm = 0.0;
	static SAMPLE *tempdatabuffer = NULL;
	int tempdatabufferlen = 0;
	if (tempdatabuffer == NULL) {
		tempdatabuffer = (SAMPLE *)malloc(data->numbytes);
	}
	while (data->isrecording) {
		;  //wait for data to stop recording
	}
	data->isprocessing = 1;  //only turn on the processing flag while copying the data.  Otherwise let the callback update the data freely.
	memcpy(tempdatabuffer, data->recordedSamples, data->numbytes);  //grab a copy of the data to process it.  During this time, recording callback will block
	tempdatabufferlen = data->maxFrameIndex;
	data->isprocessing = 0; 
	float x;

	//Ok, raw data has been collected by port audio and was just transferred to the tempdatabuffer array.  Now we get to work processing
	PlotData(&x, tempdatabuffer, tempdatabufferlen, "raw sound data");

	//First compute the 1000 sample wide moving average
	AverageFilter(tempdatabuffer, tempdatabufferlen, 1000);  //filter the data with an average filter
	PlotData(&x, tempdatabuffer, tempdatabufferlen, "filtered sound data", 2);

	//Now compute the absolute maximum of the data in the array
	int max_index;
	float max = ComputeAbsMaxOfArray(tempdatabuffer, tempdatabufferlen, &max_index);

    //Now, find all the indices (x values) of anything above 1/4 of the maximum.  This value is easily adjustable here
	std::vector<int> neighbors = FindIndicesAboveValue(tempdatabuffer, tempdatabufferlen, (SAMPLE) max*0.25);
	PlotData(neighbors, "neighbors", 3);

	//Now, find the indices of the local maximums for each group of neighbors
	std::vector<int> maxindices = FindMaxIndices(neighbors, tempdatabuffer, tempdatabufferlen);
	PlotData(maxindices, "max indices", 4);

	//Now compute the distances between the maximums
	std::vector<int> distances = FindDistances(maxindices);
	PlotData(distances, "distances", 5);

	//Now we need to do a cleanup operation on the distances to get better estimate of BPM
	distancebetweenpeaks = (float) FindMedianDistance(distances);
	std::vector<int> cleandistances = CleanUpDistances(distances, distancebetweenpeaks);
	PlotData(cleandistances, "clean distances", 6);

	//From our clean distances, compute the simple average
	distancebetweenpeaks = FindAverageDistance(cleandistances);

	//Convert to BPM via a simple formula
	bpm = ((SAMPLE) 60.0)*((SAMPLE)SAMPLE_RATE)/ distancebetweenpeaks;

	PlotEverything(tempdatabuffer, tempdatabufferlen, neighbors, maxindices, distances, cleandistances, bpm);

#if 0
	static int filesuffix = 0;
	filesuffix++;
	char filename[256];
	sprintf(filename, "audio.%d.%d", filesuffix, (int)bpm);
	FILE* fp = fopen(filename, "wb");
	if (fp) {
		int j;
		for (j = 0; j < data->maxFrameIndex; j++) {
			fprintf(fp, "%5.4f ", data->recordedSamples[j]);
		}
		fclose(fp);
	}
#endif
	return bpm;
}

/*******************************************************************/
int main(void);
int main(void)
{
	PaStreamParameters  inputParameters,
		outputParameters;
	PaStream*           stream;
	PaError             err = paNoError;
	paTestData          data;
	int                 i;
	int                 totalFrames;
	int                 numSamples;
	int                 numBytes;
	SAMPLE              max, val;
	double              average;

	printf("patest_record.c\n"); fflush(stdout);

	data.maxFrameIndex = totalFrames = NUM_SECONDS * SAMPLE_RATE; /* Record for a few seconds. */
	data.frameIndex = 0;
	numSamples = totalFrames * NUM_CHANNELS;
	numBytes = numSamples * sizeof(SAMPLE);
	data.numbytes = numBytes;
	data.isprocessing = 0;
	data.isrecording = 0;
	data.recordedSamples = (SAMPLE *)malloc(numBytes); /* From now on, recordedSamples is initialised. */
	if (data.recordedSamples == NULL)
	{
		printf("Could not allocate record array.\n");
		goto done;
	}
	for (i = 0; i < numSamples; i++) data.recordedSamples[i] = 0;

	err = Pa_Initialize();
	if (err != paNoError) goto done;

	int numDevices;
	numDevices = Pa_GetDeviceCount();
	if (numDevices < 0)
	{
		printf("ERROR: Pa_CountDevices returned 0x%x\n", numDevices);
		err = numDevices;
	}

	const   PaDeviceInfo *deviceInfo;
	for (i = 0; i < numDevices; i++)
	{
		deviceInfo = Pa_GetDeviceInfo(i);
		printf("%d %s\n", i, deviceInfo->name);
	}

	inputParameters.device = Pa_GetDefaultInputDevice(); /* default input device */
	if (inputParameters.device == paNoDevice) {
		fprintf(stderr, "Error: No default input device.\n");
		goto done;
	}
	inputParameters.channelCount = 1;                    /* mono or stereo input */
	inputParameters.sampleFormat = PA_SAMPLE_TYPE;
	inputParameters.suggestedLatency = Pa_GetDeviceInfo(inputParameters.device)->defaultLowInputLatency;
	inputParameters.hostApiSpecificStreamInfo = NULL;

	/* Record some audio. -------------------------------------------- */
	err = Pa_OpenStream(
		&stream,
		&inputParameters,
		NULL,                  /* &outputParameters, */
		SAMPLE_RATE,
		FRAMES_PER_BUFFER,
		paClipOff,      /* we won't output out of range samples so don't bother clipping them */
		recordCallback,
		&data);
	if (err != paNoError) goto done;

	err = Pa_StartStream(stream);
	if (err != paNoError) goto done;
	printf("\n=== Now recording!! Please speak into the microphone. ===\n"); fflush(stdout);

	while ((err = Pa_IsStreamActive(stream)) == 1)
	{
		Pa_Sleep(100);
		//printf("index = %d\n", data.frameIndex); fflush(stdout);
		float bpm = ComputeBPM(&data);
		printf("Current BPM: %5.4f\n", bpm);  fflush(stdout);
	}
	if (err < 0) goto done;

	err = Pa_CloseStream(stream);
	if (err != paNoError) goto done;

	/* Measure maximum peak amplitude. */
	max = 0;
	average = 0.0;
	for (i = 0; i < numSamples; i++)
	{
		val = data.recordedSamples[i];
		if (val < 0) val = -val; /* ABS */
		if (val > max)
		{
			max = val;
		}
		average += val;
	}

	average = average / (double)numSamples;

	//printf("sample max amplitude = "PRINTF_S_FORMAT"\n", max);
	printf("sample average = %lf\n", average);

	/* Write recorded data to a file. */
#if WRITE_TO_FILE
	{
		FILE  *fid;
		fid = fopen("recorded.raw", "wb");
		if (fid == NULL)
		{
			printf("Could not open file.");
		}
		else
		{
			fwrite(data.recordedSamples, NUM_CHANNELS * sizeof(SAMPLE), totalFrames, fid);
			fclose(fid);
			printf("Wrote data to 'recorded.raw'\n");
		}
	}
#endif
#if PLAYBACK
	/* Playback recorded data.  -------------------------------------------- */
	data.frameIndex = 0;

	outputParameters.device = Pa_GetDefaultOutputDevice(); /* default output device */
	if (outputParameters.device == paNoDevice) {
		fprintf(stderr, "Error: No default output device.\n");
		goto done;
	}
	outputParameters.channelCount = 2;                     /* stereo output */
	outputParameters.sampleFormat = PA_SAMPLE_TYPE;
	outputParameters.suggestedLatency = Pa_GetDeviceInfo(outputParameters.device)->defaultLowOutputLatency;
	outputParameters.hostApiSpecificStreamInfo = NULL;

	printf("\n=== Now playing back. ===\n"); fflush(stdout);
	err = Pa_OpenStream(
		&stream,
		NULL, /* no input */
		&outputParameters,
		SAMPLE_RATE,
		FRAMES_PER_BUFFER,
		paClipOff,      /* we won't output out of range samples so don't bother clipping them */
		playCallback,
		&data);
	if (err != paNoError) goto done;

	if (stream)
	{
		err = Pa_StartStream(stream);
		if (err != paNoError) goto done;

		printf("Waiting for playback to finish.\n"); fflush(stdout);

		while ((err = Pa_IsStreamActive(stream)) == 1) Pa_Sleep(100);
		if (err < 0) goto done;

		err = Pa_CloseStream(stream);
		if (err != paNoError) goto done;

		printf("Done.\n"); fflush(stdout);
	}
#endif //playback

done:
	Pa_Terminate();
	if (data.recordedSamples)       /* Sure it is NULL or valid. */
		free(data.recordedSamples);
	if (err != paNoError)
	{
		fprintf(stderr, "An error occurred while using the portaudio stream\n");
		fprintf(stderr, "Error number: %d\n", err);
		fprintf(stderr, "Error message: %s\n", Pa_GetErrorText(err));
		err = 1;          /* Always return 0 or 1, but no other return codes. */
	}
	return err;
}
