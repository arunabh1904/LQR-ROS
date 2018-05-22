#ifndef MYTIME_H 
#define MYTIME_H 


#include <string>
#include <ostream>
#include <iostream>

#define MINUTES_TO_SECONDS 60L
#define HOURS_TO_SECONDS 3600L
#define DAYS_TO_SECONDS 86400L
#define YEARS_TO_SECONDS 31536000L

using namespace std;

namespace Timing {

struct calandarDate{
		long year;
		long month;
		long day;
		long hour;
		long minute;
		long second;
		double millisecond;
};

class precisionCalendar {
public:
	precisionCalendar(const precisionCalendar* other);
	precisionCalendar(const precisionCalendar& other);
	precisionCalendar(bool getTime = true, ostream* outputStreamPointer = &cout);
	precisionCalendar( ostream* outputStreamPointer ) { precisionCalendar(true, outputStreamPointer); }
	precisionCalendar(double milliseconds, long seconds=0, long minutes=0, long hours=0, long days=0, long months=0, long years=0, ostream* outputStreamPointer = &cout);
	precisionCalendar(double milliseconds, long seconds, long minutes, long hours, long days, long months, ostream* outputStreamPointer);
	precisionCalendar(double milliseconds, long seconds, long minutes, long hours, long days, ostream* outputStreamPointer);
	precisionCalendar(double milliseconds, long seconds, long minutes, long hours, ostream* outputStreamPointer);
	precisionCalendar(double milliseconds, long seconds, long minutes, ostream* outputStreamPointer);
	precisionCalendar(double milliseconds, long seconds, ostream* outputStreamPointer);
	precisionCalendar(double milliseconds, ostream* outputStreamPointer);

	bool update();

	precisionCalendar timeSince() const;
	double secondsSince() const;


	string getTimeHMSM() const;
	string getDateTimeHMS() const;
	string getDateTimeYMDHMS() const;
	string getDateTimeYMDHMSM() const;
	
	double timeStamp() const;

	precisionCalendar &operator	+=	(double milliseconds);
	precisionCalendar &operator	-=	(double milliseconds);
	precisionCalendar operator	+	(double milliseconds);
	precisionCalendar operator	-	(double milliseconds);

	precisionCalendar &operator	+=	(const precisionCalendar& otherTime);
	precisionCalendar &operator	-=	(const precisionCalendar& otherTime);
	precisionCalendar operator	+	(const precisionCalendar& otherTime);
	precisionCalendar operator	-	(const precisionCalendar& otherTime);
	precisionCalendar operator	/	(const precisionCalendar& otherTime);
	precisionCalendar operator	*	(const precisionCalendar& otherTime);

	precisionCalendar operator*(double scale);
	precisionCalendar operator/(double scale);

	bool operator == (const precisionCalendar& otherTime) const;
	bool operator != (const precisionCalendar& otherTime) const;
	bool operator <  (const precisionCalendar& otherTime) const;
	bool operator <= (const precisionCalendar& otherTime) const;
	bool operator >  (const precisionCalendar& otherTime) const;
	bool operator >= (const precisionCalendar& otherTime) const;

	bool isPositive();

	void setTime(const precisionCalendar& other);

	string toString() const;
	operator string() const;
	operator calandarDate();
	operator double();
	friend class Timer;

protected:
	calandarDate _toDate() const;

	struct timestruct {
		long long seconds;
		double milliseconds;
	} _time;

	void _adjustForRollover();
	
	long _daysFromMonths(long month, long year, long days = 0)const;
	long _monthFromDays(long days,long year, long month = 1)const;

	bool updateMONOTONIC();

	ostream* _p_OutputStream;
};

std::ostream& operator<<(std::ostream& outputStream, precisionCalendar& time);
}
#endif // End headerfile defintion
