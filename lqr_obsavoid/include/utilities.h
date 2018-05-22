#ifndef _UTILITIES_H
#define _UTILITIES_H


#include <sstream>
#include <vector>
#include <iostream>

using namespace std;

template <class myType>
bool inRange( myType value , myType MIN, myType MAX ) {
    return (value >=MIN && value <=MAX);
}

template <class myType>
bool inRange_Array( myType value, myType MIN, myType MAX ) { 
    return (value >=MIN && value < MAX);
}


template<class myType>
string printVectorRow(const vector<myType>& myVector){
    stringstream output;
    typename vector<myType>::const_iterator myIterator;
    for(myIterator = myVector.begin(); myIterator != myVector.end(); myIterator++){
        output << *myIterator << ", ";
    }
    return output.str();
}

template<class myType>
string printVectorColumn(const vector<myType>& myVector){
    stringstream output;
    typename vector<myType>::const_iterator myIterator;
    for(myIterator = myVector.begin(); myIterator != myVector.end(); myIterator++){
        output << *myIterator << endl;
    }
    return output.str();
}

template<class myType>
string printVectorMatrix(const vector< vector< myType > >& myVector)
{
    stringstream output;
    typename vector< vector <myType> >::const_iterator outerIterators;
    vector< unsigned int > innerPosition(myVector.size(),0);
    vector< unsigned int >::iterator itP;

    bool moreToPrint = true;

    while( moreToPrint )
    {
        moreToPrint = false;

        itP = innerPosition.begin();

        for(outerIterators = myVector.begin(); outerIterators != myVector.end(); outerIterators++)
        {
            if( (*itP) < (*outerIterators).size()  )
            {
                cout << (*outerIterators)[(*itP)] << "\t";
                (*itP) ++;
                moreToPrint = true;
            } else
            {
                cout << "\t\t" ;
            }
            itP++;
        }
        cout << endl;
    }

    return output.str();
}

template<class myType>
string to_string(const myType value)
{
    stringstream output;
    output << value;
    return output.str();
}

template<typename T>
void swapPointers(T*& p_one, T* &p_two)
{
    T* tmp = p_one;
    p_one = p_two;
    p_two = tmp;
    return;
}

#endif
