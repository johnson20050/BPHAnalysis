#ifndef __histoMAP__
#define __histoMAP__

#include <map>
#include <string>
#include "TH1D.h"

class histoMAP
{
public:
    void createHisto( const std::string& name, int nbin, double min, double max );
    void copyHisto( const std::map<std::string, TH1D*>& hMap );
    void fillHisto( const std::string& name, double val );
    void clearHisto();
    void storeHistoInto( std::map<std::string, TH1D*>& totMap );
    void preventDeletedByTFile();
    ~histoMAP();

private:
    std::map<std::string, TH1D*> histoMap;
};

void histoMAP::createHisto( const std::string& name, int nbin, double min, double max )
{
    histoMap[name] = new TH1D( name.c_str(), name.c_str(), nbin, min, max );
    //histoMap[name] = new TH1D( "", "", nbin, min, max );
    return;
}
void histoMAP::copyHisto( const std::map<std::string, TH1D*>& hMap )
{
    std::map<std::string, TH1D*>::const_iterator iter = hMap.begin();    
    std::map<std::string, TH1D*>::const_iterator iend = hMap.end  ();
    while ( iter != iend )
    {
        std::pair<std::string, TH1D*> myPair( iter->first, (TH1D*)iter->second->Clone() );
        histoMap.insert(myPair);
        ++iter;
    }
    clearHisto();
    return;
}

void histoMAP::storeHistoInto( std::map<std::string, TH1D*>& totMap )
{
    std::map<std::string, TH1D*>::const_iterator tmpIt = histoMap.begin();
    std::map<std::string, TH1D*>::const_iterator tmpEd = histoMap.end  ();
    std::map<std::string, TH1D*>::      iterator totIt =   totMap.begin();
    std::map<std::string, TH1D*>::      iterator totEd =   totMap.end  ();

    // sum up the histograms to final histograms
    while ( totIt != totEd )
    {
        if ( (tmpIt = histoMap.find(totIt->first) ) != tmpEd )
            totIt->second->Add( tmpIt->second );
        ++totIt;
    }
    return;
}
void histoMAP::fillHisto( const std::string& name, double val )
{
    std::map<std::string,TH1D*>::iterator iter = histoMap.find( name );
    std::map<std::string,TH1D*>::iterator iend = histoMap.end();
    if ( iter == iend ) return;
    iter->second->Fill( val );
    return;
}
void histoMAP::clearHisto()
{
    std::map<std::string,TH1D*>::iterator iter = histoMap.begin();
    std::map<std::string,TH1D*>::iterator iend = histoMap.end()  ;
    while ( iter != iend )
        iter++->second->Reset();
    return;
}
void histoMAP::preventDeletedByTFile()
{
    std::map<std::string, TH1D*>::const_iterator iter = histoMap.begin();
    std::map<std::string, TH1D*>::const_iterator iend = histoMap.end  ();
    while ( iter != iend )
    {
        iter->second->SetDirectory(0);  // prevent histogram to be deleted by TFile::Close().
        ++iter;
    }
}
histoMAP::~histoMAP()
{
    std::map<std::string,TH1D*>::iterator iter = histoMap.begin();
    std::map<std::string,TH1D*>::iterator iend = histoMap.end()  ;
    while ( iter != iend )
        delete iter++->second;
    return;
}


void createHisto( std::map<std::string, TH1D*>& inHisto, const std::string& name, int nbin, double min, double max)
{
    inHisto[ name ] = new TH1D( name.c_str(), name.c_str(), nbin, min, max );
    //inHisto[ name ] = new TH1D( "", "", nbin, min, max );
    return;
}
void deliverHisto( const std::map<std::string, TH1D*>& tmpMap, std::map<std::string, TH1D*>& totMap )
{
    std::map<std::string, TH1D*>::const_iterator tmpIt = tmpMap.begin();
    std::map<std::string, TH1D*>::const_iterator tmpEd = tmpMap.end  ();
    std::map<std::string, TH1D*>::      iterator totIt                 ;
    std::map<std::string, TH1D*>::      iterator totEd = totMap.end  ();

    // sum up the histograms to final histograms
    while ( totIt != totEd )
    {
        if ( (tmpIt = tmpMap.find(totIt->first) ) != tmpEd )
            totIt->second->Add( tmpIt->second );
        ++totIt;
    }
    return;
}
void clearHistoMap( std::map<std::string, TH1D*>& inHistos )
{
    std::map<std::string, TH1D*>::iterator iter = inHistos.begin();
    std::map<std::string, TH1D*>::iterator iend = inHistos.end  ();
    while ( iter != iend )
        iter++->second->Reset();
    return;
}
void deleteHistoMap(std::map<std::string, TH1D*>& inHistos )
{
    std::map<std::string, TH1D*>::iterator iter = inHistos.begin();
    std::map<std::string, TH1D*>::iterator iend = inHistos.end  ();
    while ( iter != iend )
        delete iter++->second;
    return;
}
#endif
