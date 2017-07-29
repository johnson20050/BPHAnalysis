#include <iostream>
#include <stdio.h>


void eventBrowse()
{
    TFile* f1 = new TFile("/home/ltsai/ReceivedFile/DATA/testFile.root");
    TTree* t1 = new TTree("Events");
    t1->GetEntries();
}
