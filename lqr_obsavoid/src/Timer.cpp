#ifdef WIN32
#include <Windows.h>
#include <WinBase.h>
#include <MMSystem.h> // make task switching faster
#pragma comment(lib, "winmm.lib") // add library for MMSystem
#include "pthread.h"
#pragma comment(lib, "pthreadVC2.lib")
#pragma comment(lib, "pthreadVCE2.lib")
#pragma comment(lib, "pthreadVSE2.lib")
#endif

#ifdef __GNUC__
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <errno.h>
#include <string.h> // includes strerror()
#include <iostream>
#include <pthread.h>
#endif

#include<iostream>
#include <cmath>

#include "Timer.h"

using namespace std;
using namespace Timing;

void stopWatch::start(){
    startTime.update(); // resets the timer to zero
}
double stopWatch::lap(){// returns the number of seconds that have passed since start() was called;
    return startTime.secondsSince();
}
double stopWatch::stop(){// returns the number of seconds that have passed since start() was called and resets timer to zero
    double timeSince = startTime.secondsSince();
    startTime.update();
    return timeSince;
}
stopWatch::stopWatch(){
    startTime = precisionCalendar();
}
stopWatch::stopWatch(const stopWatch&  other){
    startTime  = other.startTime;
}

Timer::Timer(ostream* outputStreamPointer ):precisionCalendar(outputStreamPointer){
#ifdef WIN32
    timecaps_tag timeInfo;
    timeGetDevCaps(&timeInfo,8);
    timeBeginPeriod(timeInfo.wPeriodMin); // make time splicing as fast as possible for finer timing control
#endif
    this->updateMONOTONIC();// get realtime clock
    ticking = false;
    return;
}
Timer::Timer(double seconds, double milliseconds, ostream* outputStreamPointer):precisionCalendar(milliseconds,seconds,outputStreamPointer){
#ifdef WIN32
    timecaps_tag timeInfo;
    timeGetDevCaps(&timeInfo,8);
    timeBeginPeriod(timeInfo.wPeriodMin); // make time splicing as fast as possible for finer timing control
#endif

    ticking = false;
    return;
}
Timer::~Timer() {
    if (ticking)
        stopTicking();

#ifdef WIN32
    timecaps_tag timeInfo;
    timeGetDevCaps(&timeInfo,8);
    timeEndPeriod(timeInfo.wPeriodMax); // relqish timing request
#endif
    return;
}

bool Timer::getTime(){
    return updateMONOTONIC();
}// stores the current time in the Timer class

/**
 * Updates the stored time to the current time and returns the time elapsed.
 *
 */
Timer Timer::updateTime() {
    double timeSince = secondsSince();
    std::cout << "tSince: " << timeSince << endl;
    updateMONOTONIC();
    return Timer(timeSince,0,this->_p_OutputStream);
}

/**
 * Returns a timer that stores the elapsed time since the time that is stored in this
 * instance of the timer.
 *
 */
Timer Timer::elapsedTime() const
{

    std::cout << "Ellap Time: " << secondsSince() << std::endl;
    return Timer(0,secondsSince()/1000.0,this->_p_OutputStream);
}

void Timer::setTime(int seconds, double milliseconds) {
    _time.seconds = seconds;
    _time.milliseconds = milliseconds;
    _adjustForRollover();
}
bool Timer::waitUntilPassed()
{
    return waitTill(*this);
}
/**
 * Holds the calling thread until the given milliseconds has passed.
 *
 */
bool Timer::waitUntilPassed(double milliseconds) {
    increment(milliseconds);
    return waitUntilPassed();
}

void Timer::increment(double milliseconds){ // adds time to the stored time
    _time.milliseconds += milliseconds;
    _adjustForRollover();
}
void Timer::decrement(double milliseconds){ // subtracts time from the stored time
    _time.milliseconds -= milliseconds;
    _adjustForRollover();
}

/**
 * Stops repeatedly calling the function provided to startTicking().
 *
 */
bool Timer::stopTicking() {
    int ret = 0;

    if (!ticking)
        return true;

    ticking = false;

    if ((ret = pthread_cancel(thread_handle)) != 0) {
        std::cout << "TIMER: couldn't cancel the thread: " << strerror(ret) << ". ERROR" << std::endl;
        return false;
    }

    if ((ret = pthread_join(thread_handle, NULL)) != 0) {
        std::cout << "TIMER: there was an error stopping the thread for the timer: " << strerror(ret) << ". ERROR" << std::endl;
        return false;
    }

    return true;
}

/**
 * Returns true if the function provided to startTicking() is still being called repeatedly at the desired
 * time interval.
 *
 */
bool Timer::isTicking() const {
    return ticking;
}
/**
 * Returns the total number of seconds stored in the timer. May truncate milliseconds if seconds is too big
 *
 */
double Timer::seconds() const{
    return timeStamp();
}
/**
 * Returns the total number of seconds stored in the timer.  (Warning: you could have problems with
 * overflow if the number of seconds is too big.)
 *
 */
double Timer::milliseconds() const{
    return _time.seconds*1000.0 + _time.milliseconds;
}

bool  Timer::wait(double milliseconds){// wait for number of milliseconds

#ifdef WIN32
    long sleepTime = max(1,(long) milliseconTiming::Timer::updateTimeds);
    Sleep(sleepTime);
    return true;
#elif defined(__GNUC__)
    int ret = 0;

    timespec increment_time;
    if( clock_gettime(CLOCK_MONOTONIC, &increment_time) < 0 ) {
        std::cout << "TIMER: an error occurred when getting the current time: " + string(strerror(errno)) + ".\n";
        return false;
    }

    milliseconds += increment_time.tv_nsec/1000000L; // add current nano seconds to requested time

    increment_time.tv_sec += (long)milliseconds/1000L; // add roll-over to seconds
    milliseconds -= ((long)(milliseconds/1000L))*1000L; // subtract roll-over seconds from milliseconds

    increment_time.tv_nsec = (milliseconds)*1000000L; // add back remaining nanoseconds

    ret = clock_nanosleep(CLOCK_MONOTONIC,TIMER_ABSTIME,&increment_time,NULL);
    while( ret == EINTR ) // Sleep was interupted by a signal handeler
        ret = clock_nanosleep(CLOCK_MONOTONIC,TIMER_ABSTIME,&increment_time,NULL);

    if (ret) {
        std::cout << "TIMER: an error occurred while waiting for time to pass: " + string(strerror(ret)) + ".\n";
        return false;
    }
#endif
    return true;
}
bool  Timer::wait(const precisionCalendar& timeToWait){ // wait for the length seconds and milliseconds specified
    return wait(timeToWait._time.seconds*1000.0 + timeToWait._time.milliseconds);
}
bool  Timer::waitTill(const precisionCalendar& otherTime){// blocks until the stored time has passed
    precisionCalendar difference(otherTime);
    difference -= precisionCalendar();
    if(difference.isPositive()) {// time to wait
        return wait(difference._time.seconds*1000.0 + difference._time.milliseconds);
    }
    return true;
}

bool  Timer::waitTill(const Timer& otherTime)
{// blocks until the stored time has passed
#ifdef _WIN32
    precisionCalendar difference(otherTime);
    difference -= precisionCalendar();
    if(difference.isPositive()) {// time to wait
        return wait(difference._time.seconds*1000.0 + difference._time.milliseconds);
    }
#else
#ifdef __GNUC__ // LinuxMONOTONIC
    int ret = 0;
    timespec increment_time = otherTime.as_timespec();
//    timespec nowTime;
//    clock_gettime(CLOCK_MONOTONIC, &nowTime);

//    std::cout << "next: " << increment_time.tv_sec << "  " << increment_time.tv_nsec << endl;
//    std::cout << "now:  " << nowTime.tv_sec << "  " << nowTime.tv_nsec << endl << endl;
    ret = clock_nanosleep(CLOCK_MONOTONIC,TIMER_ABSTIME,&increment_time,NULL);
    while( ret == EINTR ) // Sleep was interupted by a signal handeler
        ret = clock_nanosleep(CLOCK_MONOTONIC,TIMER_ABSTIME,&increment_time,NULL);

    if (ret)
    {
            std::cout << "TIMER: an error occurred while waiting for time to pass: " + string(strerror(ret)) + ".\n";
            return false;
    }

#endif
#endif
    return true;
}

#ifdef __GNUC__
timespec Timer::as_timespec() const {
    timespec time;
    time.tv_sec  = this->_time.seconds;
    time.tv_nsec = this->_time.milliseconds*1000000L;
    if(time.tv_nsec == 1000000000L) {
        time.tv_nsec = 0;
        time.tv_sec +=1;
    }
    return time;
}
#endif

Timer &Timer::operator=(double milliseconds){ // sets the stored time using milliseconds
    _time.seconds = 0;
    _time.milliseconds = milliseconds;
    _adjustForRollover();
    return *this;
}

/**
 * Begins executing the function 'signal_handler' at the given timestep. (See the comment in the Timer
 * class definition for more info.)
 *
 */
bool Timer::startTicking(bool (*signal_handler)(double), double timestep_ms) {
    if (ticking)
        return true;

    NormalFunctionThreadArguments args = {signal_handler, timestep_ms, this};

    // Allocate and copy the argument structure to the private member 'function_arguments' (the space will be deleted by runNormalFunction).
    function_arguments = new unsigned char[sizeof(NormalFunctionThreadArguments)];
    memcpy(function_arguments, (unsigned char *)&args, sizeof(NormalFunctionThreadArguments));

    int ret;
    // Start up a new thread executed which executes the runNormalFunction() function with the provided arguments.
    if ( (ret = pthread_create(&thread_handle, NULL, &Timer::runNormalFunction, (void *)function_arguments)) != 0) {
        cout << "TIMER: could not create the thread for the timer -- " << strerror(ret) << ". ERROR" << endl;
        return false;
    }

    ticking = true;

    return true;
}

/**
 * This is run by the Timer thread when the user specifies a function that is not a class member
 * to the startTicking() function.
 *
 */
void *Timer::runNormalFunction(void *args) {
    Timer start_time;

    NormalFunctionThreadArguments *my_args = (NormalFunctionThreadArguments *)args;
    double timestep_ms = my_args->timestep_ms;
    Timer *timer = my_args->timer;
    bool (*signal_handler)(double) = my_args->signal_handler;

    delete[] (unsigned char *)args;

    // Store the current time.
    start_time.getTime();
    timer->getTime();

    // Keep looping until boolean ticking variable is set to false, then stop.
    while (timer->isTicking() && (*signal_handler)(start_time.elapsedTime().milliseconds())) {
        timer->increment(timestep_ms); // increment the time to the time interval

        if (!timer->waitUntilPassed()) // wait until the time interval has passed
            break;
    }

    timer->ticking = false;

    pthread_exit(NULL);
}

/**
 * Print out the stored number of seconds and milliseconds to the provided output stream.
 *
 */
ostream &operator<<(ostream &stream, const Timer &timer) {

    stream << "(s: " << timer.seconds() << ", ms: " << timer.milliseconds() << ")";

    return stream;
}
