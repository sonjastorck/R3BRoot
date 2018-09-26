#include "FairLogger.h"

#include "TClonesArray.h"
#include "FairRootManager.h"
#include "R3BLosReader.h"
#include "R3BLosMappedData.h"
extern "C" {
#include "ext_data_client.h"
#include "ext_h101_los_tamex.h"
}
#include "TMath.h"
#define IS_NAN(x) TMath::IsNaN(x)
#define NUM_LOS_DETECTORS 1
#define NUM_LOS_CHANNELS  8
#include <iostream>

using namespace std;

R3BLosReader::R3BLosReader(EXT_STR_h101_LOS_TAMEX* data, UInt_t offset)
	: R3BReader("R3BLosReader")
	, fData(data)
	, fOffset(offset)
	  , fLogger(FairLogger::GetLogger())
	  , fArray(new TClonesArray("R3BLosMappedData"))
{	
}

R3BLosReader::~R3BLosReader()
{
}

Bool_t R3BLosReader::Init(ext_data_struct_info *a_struct_info)
{
    
	int ok;

	EXT_STR_h101_LOS_TAMEX_ITEMS_INFO(ok, *a_struct_info, fOffset,
			EXT_STR_h101_LOS_TAMEX, 0);
	if (!ok) {
		perror("ext_data_struct_info_item");
		fLogger->Error(MESSAGE_ORIGIN,
				"Failed to setup structure information.");
		return kFALSE;
	}

	// Register output array in tree
	FairRootManager::Instance()->Register("LosMapped", "Land", fArray, kTRUE);
    fArray->Clear();

	// clear struct_writer's output struct. Seems ucesb doesn't do that
	// for channels that are unknown to the current ucesb config.
	EXT_STR_h101_LOS_TAMEX_onion* data = (EXT_STR_h101_LOS_TAMEX_onion*)fData;
	for (int d=0;d<NUM_LOS_DETECTORS;d++)
	{
		data->LOS[d].VTFM=0;
		data->LOS[d].VTCM=0;
		data->LOS[d].TTFLM=0;
		data->LOS[d].TTFTM=0;
		data->LOS[d].TTCLM=0;
		data->LOS[d].TTCTM=0;
	}
	return kTRUE;
}

Bool_t R3BLosReader::Read()
{
	// Convert plain raw data to multi-dimensional array
	EXT_STR_h101_LOS_TAMEX_onion* data = (EXT_STR_h101_LOS_TAMEX_onion*)fData;

	/*
	 * For variable definition, see structure EXT_STR_h101_LOS_TAMEX_onion_t 
	 * in ext_str_h101_los_tamex.h
	 *** VFTX DATA ***
	 * VTF = Size of the array TFv contaning fine VFTX times
	 * VTFv = Array containing the actual data on the fine  times
	 * VTFM = No of channels having time
	 * VTFMI = Array of TFM size containing the channel numbers of each channel with data
	 * VTFME = Array of TFM size containing the index of the first element of the next channel in data array TFv;
	 *         TFME[i]-FTME[i-1] = number of data for channel i;
	 *
	 * VTC = Size of the array TCv contaning coarse VFTX times
	 * VTCv = Array containing the actual data on the coarse times
	 * VTCM = No of channels having coarse time
	 * VTCMI = Array of TCM size containing the channel numbers of each channel with data
	 * VTCME = Array of TCM size containing the index of the first element of the next channel in data array TCv 
	 * 
	 *** TAMEX DATA *** Added by Aleksandra, Oct. 2016
	 * TTFL = Size of the array TFLv contaning fine leading times
	 * TTFLv = Array containing the actual data on the fine leading times
	 * TTFLM = No of channels having fine leading time
	 * TTFLMI = Array of TFLM size containing the channel numbers of each channel with data
	 * TTFLME = Array of TFLM size containing the index of the first element of the next channel in data array TFLv
	 *
	 * TTCL = Size of the array TCLv contaning coarse leading times
	 * TTCLv = Array containing the actual data on the coarse leading times
	 * TTCLM = No of channels having coarse leading time
	 * TTCLMI = Array of TCLM size containing the channel numbers of each channel with data
	 * TTCLME = Array of TCLM size containing the index of the first element of the next channel in data array TCLv
	 *  
	 * The same logic is for trailing times: TTFT, TTFTv, ..., TTCTMI, TTCTME
	 */

	// loop over all detectors
	
	Bool_t fprint = false;
	
	for (int d=0;d<NUM_LOS_DETECTORS;d++)
	{		    
     
     Int_t Sum = data->LOS[d].VTF+data->LOS[d].TTFT+data->LOS[d].TTFL;
     //if(data->LOS[d].TTFT != data->LOS[d].TTFL) fprint = true;
    // if(fNEvents == 9698 || fNEvents == 9701 || fNEvents == 9704) fprint = true;
		// First, we prepare time arrays for VFTX

  // VFTX first:
		uint32_t numChannels = data->LOS[d].VTFM; // not necessarly number of hits! (b/c multi hit)				
		// loop over channels
		uint32_t curChannelStart=0;     // index in v for first item of current channel
        Double_t mean_coarse_vftx = 0.;
		int sum_coarse_vftx = 0;
 // First get the average coarse time to shift all coarse counters in the same cycle (by calculateing, in the second step, deviations from the mean value.
 // If the coarse time is smaller than mean value by more than 200, then coarse counter was reseted, and thus, to its value 8192 (in case of VFTX) will be added.		
		for (int i=0;i<numChannels;i++) 
		{
			uint32_t channel = data->LOS[d].VTFMI[i]; // = 1..8
			uint32_t nextChannelStart = data->LOS[d].VTFME[i];  // index in v for first item of next channel
				
			for (int j = curChannelStart; j < nextChannelStart; j++){		
				
			int coarse_vftx = data->LOS[d].VTCv[j];
			
            mean_coarse_vftx = mean_coarse_vftx + coarse_vftx;
            sum_coarse_vftx = sum_coarse_vftx + 1;	
						
			}
			curChannelStart = nextChannelStart;
		}
            mean_coarse_vftx = mean_coarse_vftx / float(sum_coarse_vftx);


        curChannelStart=0; 
		for (int i=0;i<numChannels;i++) // VFTX, now do the mapping
		{
			uint32_t channel = data->LOS[d].VTFMI[i]; // = 1..8
			uint32_t nextChannelStart = data->LOS[d].VTFME[i];  // index in v for first item of next channel
				
			for (int j = curChannelStart; j < nextChannelStart; j++){		
				
				int coarse_vftx = data->LOS[d].VTCv[j];
			    if((mean_coarse_vftx - float(coarse_vftx)) > 200.) coarse_vftx = coarse_vftx + 8192;
			
				if(fprint) cout<<"LOS READER VFTX: "<<fNEvents<<", "<<Sum<<", "<<channel<<", "<<data->LOS[d].VTFv[j]<<", "<<data->LOS[d].VTCv[j]
				               <<", "<<coarse_vftx<<", "<<mean_coarse_vftx<<endl;					
				
				new ((*fArray)[fArray->GetEntriesFast()])				
				R3BLosMappedData(
						d+1,                  // detector number
						channel,              // channel number: 1-8
						0,                    // VFTX (0),TAMEX leading (1), TAMEX trailing (2)
						data->LOS[d].VTFv[j], // VFTX fine time
						coarse_vftx  // VFTX coarse time
						);
						
			}
			curChannelStart = nextChannelStart;
		}

// Next, TAMEX leading; first, first get the average coarse time to shift all coarse counters in the same cycle (the same as for VFTX, only here on adds 2048).     
//   if(data->LOS[d].TTFT == data->LOS[d].TTFL && data->LOS[d].TTCT == data->LOS[d].TTCL && data->LOS[d].TTFL == data->LOS[d].TTCL)
    if( 1 == 1)
	{
		numChannels = data->LOS[d].TTFLM;
		curChannelStart=0;
		Double_t mean_coarse_leading = 0.;
		int sum_coarse_leading = 0;
		for (int i=0;i<numChannels;i++)
		{
			uint32_t channel = data->LOS[d].TTFLMI[i];
			uint32_t nextChannelStart = data->LOS[d].TTFLME[i];
            
			for (int j = curChannelStart; j < nextChannelStart; j++){
				
			int coarse_leading = data->LOS[d].TTCLv[j];
			
            mean_coarse_leading = mean_coarse_leading + coarse_leading;
            sum_coarse_leading = sum_coarse_leading + 1;
            		
			}
			 
			curChannelStart = nextChannelStart;
		}  

		
		mean_coarse_leading = mean_coarse_leading / float(sum_coarse_leading);
		
		
				// Next, TAMEX leading into mapped array     
		curChannelStart=0;
		Double_t mean_coarse_lead = 0.;
		int sum_coarse_lead = 0;		
		for (int i=0;i<numChannels;i++)
		{
			uint32_t channel = data->LOS[d].TTFLMI[i];
			uint32_t nextChannelStart = data->LOS[d].TTFLME[i];
            
			for (int j = curChannelStart; j < nextChannelStart; j++){
				
			int coarse_leading = data->LOS[d].TTCLv[j];
			if((mean_coarse_leading - float(coarse_leading)) > 200.) coarse_leading = coarse_leading + 2048;
		// We now calculate again meanv alue of the "shifted" coarse leading times; this will be needed at the next step in order to
		// shift coarse trailing times in the same clock cycle as coarse leading	
            mean_coarse_lead = mean_coarse_lead + coarse_leading;
            sum_coarse_lead = sum_coarse_lead + 1;
                        	
			if(fprint) cout<<"LOS READER leading edges: "<<fNEvents<<", "<<Sum<<", "<<channel<<", "<<data->LOS[d].TTFLv[j]<<", "<<data->LOS[d].TTCLv[j]<<
			                 ", "<<coarse_leading<<"; "<<mean_coarse_leading<<endl;
				
				new ((*fArray)[fArray->GetEntriesFast()])				
				R3BLosMappedData(
						d+1,
						channel,
						1,
						data->LOS[d].TTFLv[j],
						coarse_leading
						);
		
			}
			 
			curChannelStart = nextChannelStart;
		}  
		mean_coarse_lead = mean_coarse_lead / float(sum_coarse_lead);

		// At last, TAMEX trailing: for each mapped leading edge, look for correspondiing trailing edge
	
		numChannels = data->LOS[d].TTFTM;
				 
		curChannelStart=0;
		
		for (int i=0;i<data->LOS[d].TTFTM;i++) 
		{
			
			uint32_t channel = data->LOS[d].TTFTMI[i];
			uint32_t nextChannelStart = data->LOS[d].TTFTME[i];
		
			for (int j = curChannelStart; j < nextChannelStart; j++){
			
			int coarse=data->LOS[d].TTCTv[j];
			
			if(coarse <= 25 && mean_coarse_lead >= 2023) coarse = coarse + 2048;		
			 	
			if(fprint) cout<<"LOS Reader trailing before sorting: "<<fNEvents<<", "<<Sum<<", "<<channel<<", "<<data->LOS[d].TTFTv[j]<<"; "<<coarse<<", "<<mean_coarse_lead <<endl;	
						 
			 R3BLosMappedData* mapped=NULL;
			 		 
			 Int_t n=fArray->GetEntriesFast();
			 for (int k=0;k<n;k++)
			 {
		  // Get the leading coarse time from mapped array:	
			 R3BLosMappedData* hit = (R3BLosMappedData*)fArray->At(k);
			
			 UInt_t iTypeL = hit->GetType();
			 UInt_t iCha  = hit->GetChannel();
			 
			 
			 if(iTypeL == 1) //(iCha == channel)
			 {
			 int coarse_leading=0;
			 if(iTypeL == 1) coarse_leading = hit->fTimeCoarse;	
			 
			 int tot=coarse - coarse_leading;
			 
			// if(fprint) cout<<"In trailing part: "<<n<<", "<<k<<", "<<channel<<", "<<iCha<<","<<iTypeL<<", "<<coarse_leading<<", "<<coarse<<", "<<tot<<endl;		
			 
			  			
			 if ((tot<=25) && (tot>=0))    // there is a corresponding leading edge
			  {
				//mapped=hit;
				
				new ((*fArray)[fArray->GetEntriesFast()])				
				R3BLosMappedData(
						d+1,
						channel,
						2,
						data->LOS[d].TTFTv[j],
						coarse //data->LOS[d].TTCTv[j]
						);
						
						if(fprint) cout<<"LOS Reader trailing edges: "<<fNEvents<<", "<<Sum<<", "<<channel<<", "<<data->LOS[d].TTFTv[j]<<"; "<<coarse<<", "<<tot<<endl;
								
					break;
			  }
		  }
		     }
							
			}
			curChannelStart = nextChannelStart;
		}     
	}
	else
	{
		
	}


/*		numChannels = data->LOS[d].TTFTM;
				 
		curChannelStart=0;
		
		for (int i=0;i<data->LOS[d].TTFTM;i++) 
		{
			
			uint32_t channel = data->LOS[d].TTFTMI[i];
			uint32_t nextChannelStart = data->LOS[d].TTFTME[i];
		
			for (int j = curChannelStart; j < nextChannelStart; j++){
			
			int coarse=data->LOS[d].TTCTv[j];
			
			//if(coarse == 0 && (data->LOS[d].TTCTv[j-1] >= 2047  || data->LOS[d].TTCTv[j+1] >= 2047)) coarse = 2048;
			
			//if(coarse == 1 && (data->LOS[d].TTCTv[j-1] >= 2047  || data->LOS[d].TTCTv[j+1] >= 2047)) coarse = 2049;
			
			 	
			if(fprint) cout<<"LOS Reader trailing before sorting: "<<fNEvents<<", "<<Sum<<", "<<channel<<", "<<data->LOS[d].TTFTv[j]<<"; "<<coarse<<", "<<endl;	
						 
			 R3BLosMappedData* mapped=NULL;
			 		 
			 Int_t n=fArray->GetEntriesFast();
			 for (int k=0;k<n;k++)
			 {
		  // Get the leading coarse time from mapped array:	
			 R3BLosMappedData* hit = (R3BLosMappedData*)fArray->At(k);
			
			 UInt_t iTypeL = hit->GetType();
			 UInt_t iCha  = hit->GetChannel();
			 int coarse_leading=0;
			 if(iTypeL == 1 && iCha == channel) coarse_leading = hit->fTimeCoarse;	
			 
			 //cout<<"In trailing part: "<<n<<", "<<k<<", "<<", "<<iCha<<","<<coarse_leading<<endl;		
			 int tot=coarse - coarse_leading;
							
			  if ((tot<=25) && (tot>0))    // there is a corresponding leading edge
			  {
				//mapped=hit;
				
				new ((*fArray)[fArray->GetEntriesFast()])				
				R3BLosMappedData(
						d+1,
						channel,
						2,
						data->LOS[d].TTFTv[j],
						coarse //data->LOS[d].TTCTv[j]
						);
						
						if(fprint) cout<<"LOS Reader trailing edges: "<<fNEvents<<", "<<Sum<<", "<<channel<<", "<<data->LOS[d].TTFTv[j]<<"; "<<coarse<<", "<<tot<<endl;
								
					break;
			  }
		     }
							
			}
			curChannelStart = nextChannelStart;
		}
		*/	
		 
	}     
	
	fNEvents += 1;
	     
	return kTRUE;
}

void R3BLosReader::FinishTask()
{        	
}

void R3BLosReader::Reset()
{
    // Reset the output array
    fArray->Clear();
}

ClassImp(R3BLosReader)

