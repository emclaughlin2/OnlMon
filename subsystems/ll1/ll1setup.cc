#include "LL1HEADER.h"

#include <onlmon/OnlMon.h>

// Event library

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/msg_profile.h>

// ROOT

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

void ll1setup(Event *evt, LL1HEADER *&ll1h)
{
  for ( Int_t ixmit = 0; ixmit < NXMIT; ixmit++ ) {

    Packet* p;
    p = evt->getPacket( PACKET[ixmit] );
    if ( p ) {
      ll1h->runnumber = evt->getRunNumber();
      ll1h->evtnr = p->iValue(0,"EVTNR");
      ll1h->clock = p->iValue(0,"CLOCK");

      int indx1=0; int inhit1=0; 
      int indx2=0; int inhit2=0;
      for(int is = 0; is < NRSAM; is++){
        ll1h->nhit_s1[is] =  p->iValue(is,NHITCHANNEL+NADCSH*0);
        ll1h->nhit_s2[is] =  p->iValue(is,NHITCHANNEL+NADCSH*1);
        ll1h->nhit_n1[is] =  p->iValue(is,NHITCHANNEL+NADCSH*2);
        ll1h->nhit_n2[is] =  p->iValue(is,NHITCHANNEL+NADCSH*3);

        ll1h->nhit_n[is] = ll1h->nhit_n1[is] + ll1h->nhit_n2[is];
        ll1h->nhit_s[is] = ll1h->nhit_s1[is] + ll1h->nhit_s2[is];

        ll1h->timesum_s1[is] =  p->iValue(is,NHITCHANNEL+1+NADCSH*0);
        ll1h->timesum_s2[is] =  p->iValue(is,NHITCHANNEL+2+NADCSH*1);
        ll1h->timesum_n1[is] =  p->iValue(is,NHITCHANNEL+3+NADCSH*2);
        ll1h->timesum_n2[is] =  p->iValue(is,NHITCHANNEL+4+NADCSH*3);

        ll1h->timesum_s[is] = ll1h->timesum_n1[is] + ll1h->timesum_n2[is];
        ll1h->timesum_n[is] = ll1h->timesum_s1[is] + ll1h->timesum_s2[is];

        if( ll1h->nhit_n[is]>0){
          if(inhit1 <  ll1h->nhit_n[is]){
            indx1 = is;
            inhit1 = ll1h->nhit_n[is];
          } 
        }

        if( ll1h->nhit_s[is]>0){
          if(inhit2 <  ll1h->nhit_s[is]){
            indx2 = is;
            inhit2 = ll1h->nhit_s[is];
          } 
        }

        ll1h->chargesum_n1[is] = 0;
        ll1h->chargesum_n2[is] = 0;
        ll1h->chargesum_s1[is] = 0;
        ll1h->chargesum_s2[is] = 0;
        ll1h->chargesum_s[is] = 0;
        ll1h->chargesum_n[is] = 0;

        for(int ic = 0; ic <NCH; ic++) {
          ll1h->channel[ic][is] = p->iValue(is,ic);
          if(ic<NHITCHANNEL) ll1h->chargesum_s1[is] += p->iValue(is,ic);
          else if(ic>=NADCSH*1 && ic<NHITCHANNEL+NADCSH*1) ll1h->chargesum_s2[is] += p->iValue(is,ic);
          else if(ic>=NADCSH*2 && ic<NHITCHANNEL+NADCSH*2) ll1h->chargesum_n1[is] += p->iValue(is,ic);
          else if(ic>=NADCSH*3 && ic<NHITCHANNEL+NADCSH*3) ll1h->chargesum_n2[is] += p->iValue(is,ic);
          ll1h->chargesum_n[is] = ll1h->chargesum_n1[is] + ll1h->chargesum_n2[is];
          ll1h->chargesum_s[is] = ll1h->chargesum_s1[is] + ll1h->chargesum_s2[is];
        }
        for(int it=NCH; it<(NCH+NTRIGWORDS); it++){
          ll1h->triggerwords[it - NCH][is] =  p->iValue(is,it);
        }
      }	
      ll1h->idxhitn = indx1;
      ll1h->idxhits = indx2;

      ll1h->idxsample = (ll1h->idxhitn==ll1h->idxhits) ? ll1h->idxhitn : -1;

      delete p;
    }
  }
}
