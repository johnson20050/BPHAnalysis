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


// this class record the vtxprb and event Entry
// And it becomes a cut set
// if "fillNum" is set to 0, this selection doesn't work.
// and just return true.
namespace myDef
{ typedef double myVtxprob; }
class tmpVtxprobCollection
{
public:

    tmpVtxprobCollection(unsigned fNum):fillNum(fNum) { eventNumber = 0; _startToFill = false; }
    ~tmpVtxprobCollection() {}
    void clear();
    void fillEntry( const std::pair<myDef::myVtxprob,unsigned>& dataPair, unsigned evNum );
    bool stillInCollection(unsigned evNum);
    std::map<myDef::myVtxprob, unsigned, std::greater<myDef::myVtxprob>>* getEntry();
    bool startToFill();
    unsigned getFillNumber() { return fillNum; }


private:
    // the map is sorted by "key", from bigger key to smaller key
    std::map< myDef::myVtxprob, unsigned, std::greater<myDef::myVtxprob> > valueCollection;
    std::pair<myDef::myVtxprob, unsigned> tmpPair;
    unsigned eventNumber;
    unsigned fillNum;
    bool _startToFill;
};

void tmpVtxprobCollection::clear()
{
    valueCollection.clear();
    _startToFill=false;
    valueCollection.insert(tmpPair);
}
void tmpVtxprobCollection::fillEntry( const std::pair<myDef::myVtxprob,unsigned>& dataPair, unsigned evNum ) 
{ 
    if ( stillInCollection(evNum) )
        valueCollection.insert(dataPair); 
    else
    {
        eventNumber = evNum;
        _startToFill = true;
        tmpPair = dataPair;
    }
}
// if false, fill the valueCollection! if true, continue to collect data in the same event.
bool tmpVtxprobCollection::stillInCollection(unsigned evNum)
{ return (this->eventNumber == evNum); }
std::map<myDef::myVtxprob, unsigned, std::greater<myDef::myVtxprob>>* tmpVtxprobCollection::getEntry()
{ return &valueCollection; }
bool tmpVtxprobCollection::startToFill()
{ 
    if ( fillNum == 0 ) 
        return true;
    return _startToFill; 
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

