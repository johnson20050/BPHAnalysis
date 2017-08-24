#ifndef __EVENTNUM_FILL__
#define __EVENTNUM_FILL__

// select events by vtxprob:
// Each event may not have more than 2 Lambda_b.
// so this class enable the particles filled in the histogram are restricted by event number.

// Each event can fill in up to "fillNum" particles.
// And if there are N particles in one event, sort the particles by vtxprob
// and take "fillNum" of the most vtxprob value of the particles.
#include "TH1.h"
#include <map>
#include <algorithm>
#include <functional>


template<typename T>
class tmpMassCollection
{
public:
    bool startToFill;

    tmpMassCollection(TH1* inHisto, unsigned fNum):histo(inHisto), fillNum(fNum) { eventNumber = 0; startToFill = false; }
    ~tmpMassCollection() { histoFillValue(); histo = NULL; }
    void clear();
    void fill( const std::pair<T,T>& dataPair, unsigned evNum ) 
    { 
        if ( stillInCollection(evNum) )
            valueCollection.insert(dataPair); 
        else
        {
            eventNumber = evNum;
            startToFill = true;
            tmpPair = dataPair;
        }
    }
    bool stillInCollection(unsigned evNum);
    void histoFillValue();
    void histoFillAll();

private:
    // the map is sorted by "key", from bigger key to smaller key
    std::map< T, T, std::greater<T> > valueCollection;
    std::pair<T, T> tmpPair;
    unsigned eventNumber;
    TH1* histo;
    unsigned fillNum;
};

template<typename T>
void tmpMassCollection<T>::clear()
{
    valueCollection.clear();
    startToFill=false;
    valueCollection.insert(tmpPair);
}

template<typename T>
void tmpMassCollection<T>::histoFillValue()
{
    // fill first "fillNum" particles.
    // if the candidates number is smaller than "fillNum", all of them are filled.
    unsigned numFilled = 0;
    std::map<double, double>::const_iterator iter = valueCollection.begin();
    std::map<double, double>::const_iterator iend = valueCollection.end  ();
    while ( iter != iend )
    {
        histo->Fill( iter->second );
        iter++;
        if ( ++numFilled == fillNum ) break;
    }
    this->clear();
    return;
}
template<typename T>
void tmpMassCollection<T>::histoFillAll()
{
    std::map<double, double>::const_iterator iter = valueCollection.begin();
    std::map<double, double>::const_iterator iend = valueCollection.end  ();
    while ( iter != iend )
        histo->Fill( iter++->second );

    this->clear();
    return;
}

// if false, fill the valueCollection! if true, continue to collect data in the same event.
template<typename T>
bool tmpMassCollection<T>::stillInCollection(unsigned evNum)
{ return (this->eventNumber == evNum); }



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

