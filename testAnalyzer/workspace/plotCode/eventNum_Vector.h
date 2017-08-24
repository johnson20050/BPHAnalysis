#ifndef __EVENTNUM_FILL__
#define __EVENTNUM_FILL__
#include "TH1.h"
#include <vector>
#include <algorithm>
#include <functional>

template<typename T>
class tmpMassCollection
{
public:
    std::vector<T> valueColection;
    unsigned eventNumber;
    bool startToFill;
    TH1* histo;
    unsigned fillNum;

    tmpMassCollection(TH1* inHisto, unsigned fNum):histo(inHisto), fillNum(fNum) { eventNumber = 0; startToFill = false; }
    ~tmpMassCollection() { histoFillValue(); }
    void clear() { valueColection.clear(); startToFill=false; }
    void fill( T data, unsigned evNum ) 
    { 
        if ( stillInCollection(evNum) )
            valueColection.push_back(data); 
        else
            startToFill = true;
    }
    void histoFillValue();
    bool stillInCollection(unsigned evNum);

};


void tmpMassCollection::histoFillValue()
{
    std::sort( valueColection.begin(), valueColection.end(), std::greater<T>() ); 
    unsigned vectorSize = valueColection.size();
    if ( vectorSize < fillNumber )
    {
        printf("---Warning, candidates are less than expected number! all of them filed---\n");
        for ( unsigned i = 0; i < vectorSize; ++i )
            histo->Fill( valueColection[i] );
        this->clear();
        return;
    }
    for ( unsigned i = 0; i < fillNumber; ++i )
        histo->Fill( valueColection[i] );
    this->clear();
}

// if false, fill the valueColection! if true, continue to collect data in the same event.
bool tmpMassCollection::stillInCollection(unsigned evNum)
{
    if ( this->eventNumber == evNum ) return true;
    this->eventNumber = evNum;
    startToFill=true;
    return false;
}
    
    

// prototype functions
// sort value in vector from bigger value to smaller value.
template<typename T>
void sortVector( std::vector<T>& inVector )
{ std::sort( inVector.begin(), inVector.end(), std::greater<T>() ); }

template<typename T>
void histoFillValue(TH1* histo, const std::vector<T>& inVector, unsigned fillNumber)
{
    unsigned vectorSize = inVector.size();
    if ( vectorSize < fillNumber )
    {
        printf("---Warning, candidates are less than expected number! all of them filed---\n");
        for ( unsigned i = 0; i < vectorSize; ++i )
            histo->Fill( inVector[i] );
        return;
    }
    for ( unsigned i = 0; i < fillNumber; ++i )
        histo->Fill( inVector[i] );
}






#endif
