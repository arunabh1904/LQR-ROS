/****************************************************************************
 * Timer.h  Contains the following classes
 *
 * stopWatch  // designed to track the timing of functions (real world time not processor run time) based on the system clock
 *
 * Timer	   // desined to provide precision sleep functionality
 *
 * This class can be used for timing in two different ways: 1) it can be used to call a provided function
 * at regular time intervals, and 2) it can be used for general timing.  In the first form of usage, the
 * class becomes a thread and repeatedly calls a provided function at a given time interval.  In the second
 * form of usage, the class basically becomes a stopwatch and a structure to store time.  In this form,
 * the Timer class stores time in the form of seconds (an integer with the same name) plus milliseconds
 * (a double with the same name).  Time can be added/subtracted to/from the time stored in the class. In this
 * form there are functions to find elapsed time or wait until time has passed.  The Timer class cannot
 * be used as a stopwatch and an interval timer at the same time.
 *
 * This function uses the high resolution realtime timers provided by Linux and your object file needs to be
 * compiled with the -lrt option, it also uses pthreads so your object file also needs to include the -pthread
 * option when you compile.
 *
 * Common uses...
 *
 * 1) Timing a section of code:
 *
 *     Timer my_timer;
 *     // do code here
 *     cout << my_timer.elapsedTime().milliseconds() << endl;  // report milliseconds elapsed
 *     // or...
 *     cout << my_timer.elapsedTime().seconds() << endl; // if you prefer the seconds elapsed to be reported instead
 *
 * 2) Calling a fuction, that is NOT a member of a class, at a constant time interval:
 *
 *     bool myfunction(double elapsed_time) { // print out the elapsed time since the function started getting called
 *        cout << "Hello World " << elapsed_time << endl;
 *        return true;
 *     }
 *
 *     // somewhere else in your code...
 *     Timer my_timer;
 *     my_timer.startTicking(&myfunction, 1000); // call myfunction() every 1000 milliseconds
 *     // do other things here...
 *     my_timer.stopTicking(); // stop calling myfunction()
 *
 * 3) Calling a function, that IS a member of a class, at a constant time interval:
 *
 *     class MyClass {
 *        bool mymember(double elapsed_time) { // print out the elapsed time since the function started getting called
 *           cout << "Hello World " << elapsed_time << endl;
 *           return true;
 *        };
 *     };
 *
 *     // somewhere else in your code...
 *     MyClass instance_of_myclass;
 *     Timer my_timer;
 *     my_timer.startTicking(instance_of_myclass, &MyClass::mymember, 1000); // call mymember() every 1000 milliseconds
 *     // do other things here...
 *     my_timer.stopTicking(); // stop calling myfunction()
 *
 ****************************************************************************/

#ifndef _TIMER_H
#define _TIMER_H

#include <pthread.h> // for thredding
#include <cstring>
#ifdef WIN32 // add pthread libraries
#pragma comment(lib, "pthreadVC2.lib")
#pragma comment(lib, "pthreadVCE2.lib")
#pragma comment(lib, "pthreadVSE2.lib")
#endif

#include "precisionCalendar.h"

namespace Timing{
class stopWatch{
public:
    void start(); // resets the timer to zero
    double lap(); // returns the number of seconds that have passed since start() was called;
    double stop(); // returns the number of seconds that have passed since start() or stop() was last called and resets to zero;

    stopWatch();
    stopWatch(const stopWatch&);
private:
    precisionCalendar startTime;
};

class Timer:public precisionCalendar {
private:
    bool ticking;

    void *function_arguments;

    pthread_t thread_handle;


    // This structure stores data that is passed to the timing thread when executing a callback
    // function that is a member of a class at regular intervals.
    //
    template<class T>
    struct ClassSubmemberThreadArguments {
        T *object;
        bool (T::*signal_handler)(double elapsed_ms);
        double timestep_ms;
        Timer *timer;
    };

    // This function executes the function that was provided as input to the startTicking() function
    // with the timestep interval also specified to the startTicking() function.
    //
    template<class T> static void *runClassSubmember(void *args);

    // This structure stores data that is passed to the timing thread when an ordinary function is
    // provided as a callback.
    struct NormalFunctionThreadArguments {
        bool (*signal_handler)(double elapsed_ms);
        double timestep_ms;
        Timer *timer;
    };

    static void *runNormalFunction(void *args);

public:


    Timer(ostream* outputStreamPointer  = &cout);
    Timer(double seconds, double milliseconds, ostream* outputStreamPointer  = &cout);
    ~Timer(); // destructor

        /**
         * The following functions are to be used for executing a call-back function at regular timed intervals.
         * You can't use these functions at the same time as the regular timing functions (see below).
         *
         */

        // This function starts up a new thread that repeatedly calls the provided function until stopTicking()
        // is called, or if the provided function returns false.  The provided function must return a boolean
        // and has a void argument.  The timestep_ms argument specifies the time interval to call the function at.
        // The provided function is repeatedly called until stopTicking() is called or the provided function
        // returns false.
        //
        bool startTicking(bool (*signal_handler)(double), double timestep_ms);

        // This function starts up a new thread whose sole purpose is to repeatedly call the provided function
        // that is a member of the class 'object' every 'timestep_ms' milliseconds, until the Timer's destructor
        // is called or stopTicking() is called.  If you have an instance of a AnyClass class called 'foo' and it
        // has a member called runMe(), then the appropriate usage of this function to execute runMe() every 1000
        // milliseconds is: timerclass.startTicking(foo, &AnyClass::runMe, 1000).
        //
        template<class T>
        bool startTicking(T &object, bool (T::*signal_handler)(double), double timestep_ms);

    bool stopTicking(); // stops repeatedly calling the function provided to startTicking()
    bool isTicking() const; // checks to see if the function provided to startTicking() is still being called at the desired time interval

    /**
     * These functions are to be used for your simple timing needs.  When using these functions the Timer
     * class basically becomes a storage for time that holds seconds and milliseconds.  Time can be added
     * to the Timer by using the increment() and decrement() functions or the overloaded operators.  See
     * the functions below for the available functionality and usage.
     *
     */

    bool getTime(); // stores the current time in the Timer class
    void setTime(int seconds, double milliseconds);
    Timer elapsedTime() const;// returns a Timer class containing the time elapsed since the time stored in this instance of the Timer
    Timer updateTime(); // returns the elapsed time since the time stored in this instance of the Timer and then updates the stored time
    double seconds() const;
    double milliseconds() const;

    bool waitUntilPassed();// blocks until the stored time has passed
    bool waitUntilPassed(double milliseconds); // blocks until the given time has passed, increments the timer by 'milliseconds'.

    void increment(double milliseconds); // adds time to the stored time
    void decrement(double milliseconds); // subtracts time from the stored time

    static bool wait(double milliseconds); // wait for number of milliseconds
    static bool wait(const precisionCalendar& timeToWait); // wait for the length seconds and milliseconds specified
    static bool waitTill(const precisionCalendar& otherTime); // blocks until the stored time has passed
    bool waitTill(const Timer& otherTime);// blocks until the stored time has passed

#ifdef __GNUC__
    timespec as_timespec() const;
#endif

    Timer &operator=(double milliseconds); // sets the stored time using milliseconds

    friend ostream &operator<<(ostream &stream, const Timer &time);
};
ostream &operator<<(ostream &stream, const Timer &timer); // print out the stored time to stdout

/**
 * Begins executing the function 'signal_handler' that is a member of the 'object' class repeatedly
 * at the desired timestep.  (See the comment in the Timer class definition for more info.)
 *
 */
template<class T>
bool Timer::startTicking(T &object, bool (T::*signal_handler)(double), double timestep_ms) {
    int ret = 0;

    if (ticking)
        return true;

    // Create a new instance of the thread argument structure and populate it with information.
    ClassSubmemberThreadArguments<T> args = {&object, signal_handler, timestep_ms, this};

    // Allocate and copy the argument structure to the private member 'function_arguments' (the space will be deleted by runClassSubmember).
    function_arguments = new unsigned char[sizeof(ClassSubmemberThreadArguments<T>)];
    memcpy(function_arguments, (unsigned char *)&args, sizeof(ClassSubmemberThreadArguments<T>));

    thread_handle = 0;
    // Start up a new thread executed which executes the runClassSubmember() function with the provided arguments.
    if ((ret = pthread_create(&thread_handle, NULL, &Timer::runClassSubmember<T>, function_arguments)) != 0 && thread_handle!=0) {
        std::cout << "TIMER: there was an error starting the thread for the timer: " << strerror(ret) << " ERROR" << std::endl;
        return false;
    }

    ticking = true;

    return true;
}

/**
 * This is run by the Timer thread when the user specifies a function that is a member of some class
 * to the startTicking() function.
 *
 */
template<class T>
void *Timer::runClassSubmember(void *args) {


    ClassSubmemberThreadArguments<T> *my_args = (ClassSubmemberThreadArguments<T> *)args;
    T *object = my_args->object;
    bool (T::*signal_handler)(double) = my_args->signal_handler;
    double timestep_ms = my_args->timestep_ms;
    Timer *timer = my_args->timer;

    delete[] (unsigned char *)args; // delete the place in memory where the function arguments were stored

    // Store the current time.
    stopWatch start_time;
    timer->getTime();

    // Keep looping until boolean ticking variable is set to false, then stop.
    while (timer->isTicking() && (object->*signal_handler)(start_time.stop()*1000.0) )
    {
        timer->increment(timestep_ms); // increment the time to the time interval

        if (!timer->waitUntilPassed()) // wait until the time interval has passed
            break;
    }

    timer->ticking = false;

    pthread_exit(NULL);
}

}
#endif
