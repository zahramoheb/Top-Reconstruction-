#define MyClass_cxx
#include "MyClass.h"
#include <iostream>
#include <math.h>
using namespace std;

TH1D *hh1= new TH1D("Pz Reconstruction Error", "Pz Reconstruction Error", 100, -3, 3);
TH1D *hh2 = new TH1D("leptonic invariant mass of top quark", " leptonic invariant mass of top quark", 100,-100, 300);
TH1D *hh3 = new TH1D("hadronic invariant mass of top quark", "hadronic invariant mass of top quark", 100, 100, 300);

int main()
{
  MyClass LHEFIPM;
 
  // Event_lightquark_num & Event_lepton_num count the number of lightquark & lepton in the event
  int Event_lightquark_num, Event_lepton_num;
  
  // noulepton_pos, lepton_pos, b_pos, bbar_pos, q_pos ,and qbar_pos save the position (0-11) of particles in the event
  int noulepton_pos, lepton_pos, b_pos, bbar_pos, q_pos, qbar_pos;
  
  // numberofsemi counts total number of semileptonic Event in all events.
  int numberofsemi=0;
  
  //ax^2 + bx + c = 0, delta=b^2-4a*c, X1 and X2 are the answers of the equation, Pz_rec is reconstruction of noulepton Pz
  double_t equ_a = 1, equ_b, equ_c, delta, X1, X2, Pz_rec;
  
  //error_rec is error of Pz reconstruction
  double_t error_rec;
  
  //w_invar is constant, c_speed is constant
  double_t w_invar=80.379, c_speed= 1 ;
  
  // leptonic invariant mass of top quark
  double_t lep_topmass;
  
  // hadronic invariant mass of top quark
  double_t had_topmass;

  // create TCanvas
  TCanvas *cc1 = new TCanvas("Top Quark Reconstruction", "Top Quark Reconstruction", 1000, 800);
 
  // loop over all Events
  for(Int_t i=0; i < LHEFIPM.EvGetEntries() ; i++)
    {
    
    LHEFIPM.GetEntry(i);
   
    // reset counter
    Event_lightquark_num = 0 ;
    Event_lepton_num = 0 ;
   
    // loop over the Event Particles   
    for(Int_t j=0; j < LHEFIPM.kMaxParticle ; j++)
      {
	
      // Check stable particle
      if(LHEFIPM.Particle_Status[j] == 1)
	{
	  
	// Check PID for q lightquark
	if(LHEFIPM.Particle_PID[j] == 1 || LHEFIPM.Particle_PID[j] == 2 ||  LHEFIPM.Particle_PID[j] == 3 )
	  {
	    
	  //count lightquark
	  Event_lightquark_num++;
	  
	  //save the position of q
	  q_pos = j;
	  
	  }
	
	// Check PID for qbar lightquark
	if( LHEFIPM.Particle_PID[j] == -1 ||LHEFIPM.Particle_PID[j] == -2 || LHEFIPM.Particle_PID[j] == -3)
	  {
	    
	  //count lightquark
	  Event_lightquark_num++;
	  
	  //save the position of qbar
	  qbar_pos = j;
	  }
	
	// Check PID for lepton
	if(LHEFIPM.Particle_PID[j] == 11 || LHEFIPM.Particle_PID[j] == -11 || LHEFIPM.Particle_PID[j] == 13 || LHEFIPM.Particle_PID[j] == -13 || LHEFIPM.Particle_PID[j] == 15 || LHEFIPM.Particle_PID[j] == -15 || LHEFIPM.Particle_PID[j] == 17 || LHEFIPM.Particle_PID[j] == -17 )
	  {
	    
	  //count lepton
	  Event_lepton_num++;
	  
	  //save the position of lepton , it just works for semileptonic and for fullleptonic needs lepton_pos array
	  lepton_pos = j;
	  
	  }
	// Check PID for noulepton
	if( LHEFIPM.Particle_PID[j] == 12 || LHEFIPM.Particle_PID[j] == -12 || LHEFIPM.Particle_PID[j] == 14 || LHEFIPM.Particle_PID[j] == -14 ||  LHEFIPM.Particle_PID[j] == 16 || LHEFIPM.Particle_PID[j] == -16 || LHEFIPM.Particle_PID[j] == 18 || LHEFIPM.Particle_PID[j] == -18)
	  {
	    
	  //count lepton
	  Event_lepton_num++;
	  
	  //save the position of noulepton , it just works for semileptonic and for fullleptonic needs noulepton_pos array
	  noulepton_pos = j;
	   
          }
	
	// Check PID for b
	if( LHEFIPM.Particle_PID[j] == 5)
	  {
	    
	  //save the position of noulepton
	  b_pos=j;
	  
	  }
	
	// Check PID for bbar
	if( LHEFIPM.Particle_PID[j] == -5)
	  
	  {
	    
	  //save the position of noulepton
	  bbar_pos=j;
	  
	  }
	
	} //End of Check stable particle
      
      }//End of loop over the Event Particles
    
      // Check semileptonic event
      if(Event_lightquark_num == 2 && Event_lepton_num == 2 )
	{
	//solve equation  
	equ_b = 2 * LHEFIPM.Particle_Pz[lepton_pos];
	equ_c = -pow((LHEFIPM.Particle_E[lepton_pos])/c_speed,2)+pow((LHEFIPM.Particle_Px[lepton_pos]),2)+pow((LHEFIPM.Particle_Py[lepton_pos]),2)+pow((LHEFIPM.Particle_Pz[lepton_pos]),2)-pow((LHEFIPM.Particle_E[noulepton_pos])/c_speed , 2)+pow((LHEFIPM.Particle_Px[noulepton_pos]), 2)+pow((LHEFIPM.Particle_Py[noulepton_pos]),2)+2*(-(LHEFIPM.Particle_E[lepton_pos]/c_speed)*(LHEFIPM.Particle_E[noulepton_pos]/c_speed)+LHEFIPM.Particle_Px[lepton_pos]*LHEFIPM.Particle_Px[noulepton_pos]+LHEFIPM.Particle_Py[lepton_pos]*LHEFIPM.Particle_Py[noulepton_pos])+(w_invar)*(w_invar);
        delta = pow((equ_b),2)-(4*equ_a*equ_c);

	//check delta and choose right answer
	if(delta>0)
	  {
	    X1=(-(equ_b)-sqrt(delta))/(2*(equ_a));
	    X2=(-(equ_b)+sqrt(delta))/(2*(equ_a));
	    if(abs(X1)<abs(X2))
	      {
		Pz_rec = X1;
	      }
	    else
	      {
		Pz_rec = X2;
	      }
	  }
	else
	  {	    
	    Pz_rec = (-(equ_b))/(2*(equ_a));
	  }
	// calculate error of reconstruction
	error_rec =((LHEFIPM.Particle_Pz[noulepton_pos])-(Pz_rec))/(LHEFIPM.Particle_Pz[noulepton_pos]);
	

	//fill hh1 with error of reconstruction
	hh1->Fill(error_rec);

	// calculate leptonic top quark mass
        lep_topmass= sqrt(pow((LHEFIPM.Particle_E[lepton_pos])/c_speed +(LHEFIPM.Particle_E[noulepton_pos])/c_speed + (LHEFIPM.Particle_E[b_pos])/c_speed,2)-pow( LHEFIPM.Particle_Px[lepton_pos] + LHEFIPM.Particle_Px[noulepton_pos] + LHEFIPM.Particle_Px[b_pos] ,2) - pow (LHEFIPM.Particle_Py[lepton_pos] + LHEFIPM.Particle_Py[noulepton_pos] + LHEFIPM.Particle_Py[b_pos] ,2)- pow( LHEFIPM.Particle_Pz[lepton_pos]+ Pz_rec + LHEFIPM.Particle_Pz[b_pos] ,2) );
	
	// fill hh2 with mass of leptonic top quark
        hh2->Fill(lep_topmass);

	// calculate hadronic top quark mass
	had_topmass = sqrt(pow((LHEFIPM.Particle_E[q_pos])/c_speed + (LHEFIPM.Particle_E[qbar_pos])/c_speed + (LHEFIPM.Particle_E[bbar_pos])/c_speed,2)-pow( LHEFIPM.Particle_Px[q_pos] + LHEFIPM.Particle_Px[qbar_pos] + LHEFIPM.Particle_Px[bbar_pos] ,2) - pow (LHEFIPM.Particle_Py[q_pos] + LHEFIPM.Particle_Py[qbar_pos] + LHEFIPM.Particle_Py[bbar_pos] ,2)- pow( LHEFIPM.Particle_Pz[q_pos]+ LHEFIPM.Particle_Pz[qbar_pos] + LHEFIPM.Particle_Pz[bbar_pos] ,2) );

	// fill hh3 with mass of hadronic top quark
        hh3->Fill(had_topmass);

	//count semileptonic
	numberofsemi++;

	}// End of Check semileptonic event
      
  }// End of loop over all Events


  
  

  // hh1->Draw();

  //  hh2->Draw();
 
     hh3->Draw();

  return numberofsemi;
}// End of main

