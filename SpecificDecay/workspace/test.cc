// this struct uses the map on particle name and reco::Candidate
// and finally pair particle name and RefCountedKinematicParticle
struct recoParticleInfo
{
    recoParticleInfo( const std::string& compName, const std::string& daugName, const reco::Candidate* inptr ) :
        cName( compName ),
        dName( daugName ),
        cptr( inptr ),
        refptr( nullptr ) {}
    std::string getFullName() const
    { 
        if ( this->isCompDaughter() )
            return cName+"/"+dName; 
        return dName;
    }
    std::string getCompName() const
    { return cName; }
    std::string getDaugName() const
    { return dName; }
    int getRecoCharge() const
    { return cptr->charge(); }
    float getRecoPt() const
    { return cptr->pt(); }
    GlobalVector getRefMom() const
    { return refptr->currentState().kinematicParameters().momentum(); }
    int getRefCharge() const
    { return refptr->currentState().particleCharge(); }
    void setRefParticle( RefCountedKinematicParticle& obj )
    { refptr.swap(obj); }
    RefCountedKinematicParticle getRefParticle() const
    { return refptr; }
    bool isCompDaughter() const
    { return !cName.empty(); }
    bool checkSameCharge() const
    {
        if ( cptr && refptr ) 
        {
            if ( this->getRecoCharge() == this->getRefCharge() )
                return true;
        }
        else
            printf ("recoParticleInfo::cptr or refptr doesn't exit!\n");
        return false;
    }
    private:        
    std::string cName;
    std::string dName;
    const reco::Candidate* cptr;
    RefCountedKinematicParticle refptr;
};
