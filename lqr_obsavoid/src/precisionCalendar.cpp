#ifdef WIN32
	#include <Windows.h>
#endif

#ifdef __GNUC__
#include <sys/time.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#endif


#include "precisionCalendar.h"
#include "utilities.h"
#include <time.h>

#include <string>
#include <iostream>
using namespace std;
using namespace Timing;

#define YEAR_ZERO 1970L

precisionCalendar::precisionCalendar(const precisionCalendar* other) {
	_time = other->_time;
	_p_OutputStream = 	other->_p_OutputStream;;
}
precisionCalendar::precisionCalendar(const precisionCalendar& other) {
	_p_OutputStream = 	other._p_OutputStream;
	_time = other._time;
}
precisionCalendar::precisionCalendar(bool getTime, ostream* outputStreamPointer) {
	_p_OutputStream = 	outputStreamPointer;

	if(getTime){
		update();
	}else{
		_time.seconds = 0;
		_time.milliseconds = 0;
	}
	_adjustForRollover();
}

precisionCalendar::precisionCalendar(double milliseconds, long seconds, long minutes, long hours, long days, long months, long years, ostream* outputStreamPointer){
	_p_OutputStream = 	outputStreamPointer;

	_time.seconds		= 0;
	_time.seconds		+= (years)*YEARS_TO_SECONDS;
	_time.seconds		+= (years)/4*DAYS_TO_SECONDS;// Leap year days
	_time.seconds		-= ((years)%4==0)*DAYS_TO_SECONDS*(years!=0);// if this is a leap year remove it
	_time.seconds		+= _daysFromMonths(months-1,years)*DAYS_TO_SECONDS;
	_time.seconds		+= (days)*DAYS_TO_SECONDS;
	_time.seconds		+= hours*HOURS_TO_SECONDS;
	_time.seconds		+= minutes*MINUTES_TO_SECONDS;
	_time.seconds		+= seconds;
	_time.milliseconds	=  milliseconds;
	
	_adjustForRollover();
}

precisionCalendar::precisionCalendar(double milliseconds, long seconds, long minutes, long hours, long days, long months, ostream* outputStreamPointer) {
	precisionCalendar( milliseconds, seconds, minutes, hours, days, months, 0 , outputStreamPointer);
}
precisionCalendar::precisionCalendar(double milliseconds, long seconds, long minutes, long hours, long days, ostream* outputStreamPointer) {
	precisionCalendar( milliseconds, seconds, minutes, hours, days, 0, 0 , outputStreamPointer);
}
precisionCalendar::precisionCalendar(double milliseconds, long seconds, long minutes, long hours, ostream* outputStreamPointer) {
	precisionCalendar( milliseconds, seconds, minutes, hours, 0 , 0 , 0 , outputStreamPointer);
}
precisionCalendar::precisionCalendar(double milliseconds, long seconds, long minutes, ostream* outputStreamPointer) {
	precisionCalendar( milliseconds, seconds, minutes, 0 , 0 , 0 , 0 , outputStreamPointer);
}
precisionCalendar::precisionCalendar(double milliseconds, long seconds, ostream* outputStreamPointer) {
	precisionCalendar( milliseconds, seconds, 0 , 0 , 0 , 0 , 0 , outputStreamPointer);
}
precisionCalendar::precisionCalendar(double milliseconds, ostream* outputStreamPointer) {
	precisionCalendar( milliseconds, 0 , 0 , 0 , 0 , 0 , 0 , outputStreamPointer);
}

calandarDate precisionCalendar::_toDate() const{
	calandarDate date;

	long long secondsLeft = _time.seconds;

	date.year = (long) secondsLeft/YEARS_TO_SECONDS;

	secondsLeft -= (date.year*YEARS_TO_SECONDS) + ((date.year+YEAR_ZERO)/4-YEAR_ZERO/4)*DAYS_TO_SECONDS*(date.year!=0)// Remove leap year days
		- (((date.year+YEAR_ZERO)%4==0)*DAYS_TO_SECONDS)*(date.year!=0);//if this is a leep year put the day back

	date.year	+= YEAR_ZERO * (date.year!=0);
	if(date.year != 0){
		date.month	= _monthFromDays((long)secondsLeft/DAYS_TO_SECONDS,date.year);
		secondsLeft -= _daysFromMonths(date.month-1,date.year)*DAYS_TO_SECONDS;
	}else{
		date.month = 0;
	}
	date.day	= 	(long)secondsLeft/DAYS_TO_SECONDS;
	secondsLeft 	-= 	date.day*DAYS_TO_SECONDS;
	date.day 	+=	1*(date.month!=0 || date.year!=0); // months dont start at day zero
	date.hour	= 	(long)secondsLeft/HOURS_TO_SECONDS;
	secondsLeft 	-= 	date.hour*HOURS_TO_SECONDS;
	date.minute 	= 	(long)secondsLeft/MINUTES_TO_SECONDS;
	secondsLeft 	-= 	date.minute*MINUTES_TO_SECONDS;
	date.second 	= 	(long)secondsLeft;

	date.millisecond = _time.milliseconds;

	return date;
}

string precisionCalendar::getDateTimeHMS() const{
	std::string output="";

	calandarDate date;
	date = _toDate();

    output += to_string<long long>( date.hour);
	output += ":";
    output += to_string<long long>(  date.minute);
	output += ":";
    output += to_string<long long>(  date.second);

	return output;
}

string precisionCalendar::getDateTimeYMDHMS() const{
	std::string output;
	calandarDate date = _toDate();
#include <time.h>
    output = to_string<long long>(  date.year);
	output += "-";
    output += to_string<long long>(  date.month);
	output += "-";
    output += to_string<long long>(  date.day);
	output += " ";
    output += to_string<long long>(  date.hour);
	output += "-";
    output += to_string<long long>(  date.minute);
	output += "-";
    output += to_string<long long>(  date.second);

	return output;
}

string precisionCalendar::getDateTimeYMDHMSM() const{
	std::string output;
	output =getDateTimeYMDHMS();
	output += ".";
	long long timeMilliseconds = (long long) _time.milliseconds;
	
    output += to_string<long long>( timeMilliseconds/100);
	timeMilliseconds -= (timeMilliseconds/100)*100;
    output += to_string<long long>( timeMilliseconds/10);
	timeMilliseconds -= (timeMilliseconds/10)*10;
    output += to_string<long long>( timeMilliseconds);

	return output;
}

string precisionCalendar::getTimeHMSM() const{
	std::string output;
	output = getDateTimeHMS();
	output += ".";
	long long timeMilliseconds = (long long)_time.milliseconds;
	
    output += to_string<long long>( timeMilliseconds/100);
	timeMilliseconds -= (timeMilliseconds/100)*100;
    output += to_string<long long>( timeMilliseconds/10);
	timeMilliseconds -= (timeMilliseconds/10)*10;
    output += to_string<long long>( timeMilliseconds);

	return output;
}

bool precisionCalendar::update(){
	#ifdef WIN32
		SYSTEMTIME systime;
		GetLocalTime(&systime);
		_time.seconds		= 0;
		_time.seconds		+= ((systime.wYear)-YEAR_ZERO)*YEARS_TO_SECONDS;
		_time.seconds		+= ((systime.wYear)-YEAR_ZERO)/4*DAYS_TO_SECONDS;// Leap year days
		_time.seconds		-= (((systime.wYear)-YEAR_ZERO)%4==0)*DAYS_TO_SECONDS;// if this is a leap year remove it
		_time.seconds		+= _daysFromMonths(systime.wMonth-1,systime.wYear)*DAYS_TO_SECONDS;
		_time.seconds		+= (systime.wDay-1)*DAYS_TO_SECONDS;
		_time.seconds		+= systime.wHour*HOURS_TO_SECONDS;
		_time.seconds		+= systime.wMinute*MINUTES_TO_SECONDS;
		_time.seconds		+= systime.wSecond;
		_time.milliseconds	=  systime.wMilliseconds;
	#endif

	#ifdef __GNUC__
		timespec  _systime;
		if( clock_gettime(CLOCK_REALTIME, &_systime) < 0 ) {
			*_p_OutputStream << "TIMER: an error occurred when getting the current time: " + string(strerror(errno)) + ".\n";
			return false;
		}
		_time.seconds = _systime.tv_sec;
		_time.milliseconds = _systime.tv_nsec/1000000.0;
	#endif
	return true;
}

bool precisionCalendar::updateMONOTONIC() {
#ifdef _WIN32
	return update();
#elif defined(__GNUC__) // linux
	timespec _systime;
	if( clock_gettime(CLOCK_MONOTONIC, &_systime) < 0 ) {
		*_p_OutputStream << "TIMER: an error occurred when getting the current time: " + string(strerror(errno)) + ".\n";
		return false;
	}
	_time.seconds = _systime.tv_sec;
	_time.milliseconds = _systime.tv_nsec/1000000.0;

	return true;
#endif
	return false;

}

precisionCalendar precisionCalendar::timeSince()const{
	precisionCalendar now;
	return now - (*this);
}
double	precisionCalendar::secondsSince() const{
	long double seconds;

	precisionCalendar now;
	seconds = (long double) now._time.seconds-this->_time.seconds;
	seconds += (now._time.milliseconds-this->_time.milliseconds)/1000.0;

	return seconds;

}
long precisionCalendar::_daysFromMonths(long month,long year, long days)const{
	
	switch(month){
	case 1:
		days += 31;
		break;
	case 2:
		if (year%4 == 0)
			days += 29;
		else
			days += 28;
		break;
	case 3:
		days += 31;
		break;
	case 4:
		days += 30;
		break;
	case 5:
		days += 31;
		break;
	case 6:
		days += 30;
		break;
	case 7:
		days += 31;
		break;
	case 8:
		days += 31;
		break;
	case 9:
		days += 30;
		break;
	case 10:
		days += 31;
		break;
	case 11:
		days += 30;
		break;
	case 12:
		days += 31;
		break;
	default:
		break;
	}

	if(month <= 1 || month > 12)
		return days;

	return _daysFromMonths( month-1, year,days);
	
}

long precisionCalendar::_monthFromDays(long days,long year, long month)const{
	switch((month)){
	case 1:
		days -= 31;
		break;
	case 2:
		if (year%4 == 0)
			days -= 29;
		else
			days -= 28;
		break;
	case 3:
		days -= 31;
		break;
	case 4:
		days -= 30;
		break;
	case 5:
		days -= 31;
		break;
	case 6:
		days -= 30;
		break;
	case 7:
		days -= 31;
		break;
	case 8:
		days -= 31;
		break;
	case 9:
		days -= 30;
		break;
	case 10:
		days -= 31;
		break;
	case 11:
		days -= 30;
		break;
	case 12:
		days -= 31;
		break;
	case 0:
		break;
	}
	if(days <= 0)
		return (month);

	return _monthFromDays( days,year, month+1);
	
}

void precisionCalendar::_adjustForRollover(){
	if( _time.milliseconds >= 1000.0 || _time.milliseconds <= -1000.0){
		_time.seconds += (long long)_time.milliseconds/1000LL;
		_time.milliseconds -= (_time.milliseconds/1000L) * 1000.0;
	}
	if(_time.seconds >= 0 && _time.milliseconds < 0.0){
		_time.seconds -=1L;
		_time.milliseconds += 1000.0;
	}else if(_time.seconds < 0 && _time.milliseconds > 0.0){
		_time.seconds +=1L;
		_time.milliseconds -= 1000.0;
	}
}

precisionCalendar	&precisionCalendar::operator+=	(double milliseconds){
		_time.milliseconds+=milliseconds;
		_adjustForRollover();
		return *this;
}
precisionCalendar	&precisionCalendar::operator-=	(double milliseconds){
		_time.milliseconds-=milliseconds;
		_adjustForRollover();
		return *this;
}
precisionCalendar	precisionCalendar::operator+	(double milliseconds){
		precisionCalendar newTime(this);
		newTime._time.milliseconds	+= milliseconds;
		newTime._adjustForRollover();

		return newTime;
	}
precisionCalendar	precisionCalendar::operator-	(double milliseconds){
		precisionCalendar newTime(this);
		newTime._time.milliseconds	-= milliseconds;
		newTime._adjustForRollover();

		return newTime;
	}

precisionCalendar	precisionCalendar::operator+	(const precisionCalendar& otherTime){
	precisionCalendar newTime(false);
	newTime._time.seconds		= _time.seconds			+ otherTime._time.seconds;
	newTime._time.milliseconds	= _time.milliseconds	+ otherTime._time.milliseconds;
	newTime._adjustForRollover();

	return newTime;
}
precisionCalendar	&precisionCalendar::operator+=	(const precisionCalendar& otherTime){
	_time.seconds		+= otherTime._time.seconds;
	_time.milliseconds	+= otherTime._time.milliseconds;
	_adjustForRollover();
	return *this;
}
precisionCalendar	precisionCalendar::operator-	(const precisionCalendar& otherTime){	
	precisionCalendar newTime(false);
	newTime._time.seconds		= _time.seconds			- otherTime._time.seconds;
	newTime._time.milliseconds	= _time.milliseconds	- otherTime._time.milliseconds;
	
	newTime._adjustForRollover();

	return newTime;

}
precisionCalendar	&precisionCalendar::operator-=	(const precisionCalendar& otherTime){
	_time.seconds		-= otherTime._time.seconds;
	_time.milliseconds	-= otherTime._time.milliseconds;
	_adjustForRollover();
	return *this;
}

bool  precisionCalendar::operator ==  (const precisionCalendar	& otherTime) const{
	return (_time.seconds == otherTime._time.seconds && _time.milliseconds == otherTime._time.milliseconds);
}
bool  precisionCalendar::operator !=  (const precisionCalendar	& otherTime) const{
	return !(*this==otherTime);
}
bool  precisionCalendar::operator <  (const precisionCalendar	& otherTime)const {

	return ((_time.seconds<otherTime._time.seconds) || (_time.seconds  == otherTime._time.seconds && _time.milliseconds<otherTime._time.milliseconds));
}
bool  precisionCalendar::operator <=	(const precisionCalendar	& otherTime)const {
	return *this < otherTime || *this == otherTime;
}
bool  precisionCalendar::operator > (const precisionCalendar & otherTime) const {
	
	return ((_time.seconds>otherTime._time.seconds) || (_time.seconds  == otherTime._time.seconds && _time.milliseconds>otherTime._time.milliseconds));

}
bool  precisionCalendar::operator >=	(const precisionCalendar & otherTime)const {

	return *this > otherTime || *this == otherTime;

}

double precisionCalendar::timeStamp() const{
	long double seconds;
	seconds = (long double) _time.seconds;
	seconds += _time.milliseconds/1000.0;
	return seconds;

}
bool precisionCalendar::isPositive() { // returns true if the time held is a positive value false otherwise
	this->_adjustForRollover();
	return (_time.seconds > 0 && _time.milliseconds > 0);
}
void precisionCalendar::setTime(const precisionCalendar& other){
	_time.milliseconds = other._time.milliseconds;
	_time.seconds = other._time.seconds;
	_adjustForRollover();
}
string precisionCalendar::toString() const{
	calandarDate date = _toDate();
	string returnString = "";
	int millisecond_hundreds = (long)date.millisecond / 100L;
	int millisecond_tens = ((long)date.millisecond-millisecond_hundreds*100L)/10;
	int millisecond_ones = ((long)date.millisecond-millisecond_hundreds*100-millisecond_tens*10);
	
	if(date.year!=0){
        returnString += to_string<long long>(  date.year) + "-" + to_string<long long>( date.month) + "-" ;
	}
    returnString += to_string<long long>( date.day) + " ";
	if(date.hour < 10)
		returnString += "0";
    returnString += to_string<long long>( date.hour) + ":" ;
	
	if(date.minute < 10)
		returnString += "0";
    returnString += to_string<long long>( date.minute) + ":" ;
	
	if(date.second < 10)
		returnString += "0";
    returnString += to_string<long long>( date.second) + ".";
	
    returnString += to_string<long long>( millisecond_hundreds) + to_string<long long>( millisecond_tens) + to_string<long long>( millisecond_ones);
	
	return returnString;
}
precisionCalendar::operator string()const{
	return toString();
}
precisionCalendar::operator calandarDate(){
	return _toDate();
}
precisionCalendar::operator double(){
	return timeStamp();
}

std::ostream& operator<<(std::ostream& outputStream, precisionCalendar& time){
	outputStream << (string)time;
	return outputStream;
}
