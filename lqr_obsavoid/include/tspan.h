#ifndef TSPAN_H
#define TSPAN_H

struct T_Span
{
    double t0;
    double tf;

    T_Span()
    {
        t0=0;
        tf=0;
    }

    T_Span(double tf_)
    {
        t0=0;
        tf=tf_;
    }
    T_Span(double to_, double tf_)
    {
        t0=to_;
        tf=tf_;
    }
};

#endif // TSPAN_H

