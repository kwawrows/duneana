////////////////////////////////////////////////////////////////////////
// Class:       TriggerPrimitiveFinderPass2
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "duneana/DAQSimAna/TriggerPrimitiveFinderPass1.h"

#include "duneana/DAQSimAna/AlgParts.h"

#include <algorithm> // for std::transform

class TriggerPrimitiveFinderPass2 : public TriggerPrimitiveFinderPass1 {
public:
  explicit TriggerPrimitiveFinderPass2(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

    virtual std::vector<TriggerPrimitiveFinderTool::Hit>
    findHits(const std::vector<unsigned int>& channel_numbers, 
             const std::vector<std::vector<short>>& collection_samples);
    

protected:

    void hitFinding(const std::vector<short>& waveform,
                    const std::vector<short>& iqr,
                    std::vector<TriggerPrimitiveFinderTool::Hit>& hits,
                    int channel);
private:
    float m_sigmaThreshold;
};


TriggerPrimitiveFinderPass2::TriggerPrimitiveFinderPass2(fhicl::ParameterSet const & p)
    : TriggerPrimitiveFinderPass1(p),
      m_sigmaThreshold(p.get<float>("ThresholdInSigma", 500))
// Initialize member data here.
{
  std::cout << "Threshold in sigma is " << m_sigmaThreshold  << " (ignore the Threshold = 10 on previous line)\n";
}

void
TriggerPrimitiveFinderPass2::hitFinding(const std::vector<short>& waveform,
                                        const std::vector<short>& iqr,
                                        std::vector<TriggerPrimitiveFinderTool::Hit>& hits,
                                        int channel)
{
    //---------------------------------------------
    // Hit finding
    //---------------------------------------------
    bool is_hit=false;
    bool was_hit=false;
    std::vector<int> hit_charge; 

    //Initialize current hit
    TriggerPrimitiveFinderTool::Hit hit(channel, 0, 0, 0, 0, 0);

    //Loop over ADCs in the waveform 
    for(size_t isample=0; isample<waveform.size()-1; ++isample){

        int sample_time=isample*m_downsampleFactor;
        short adc=waveform[isample];
        is_hit=(float)adc>m_sigmaThreshold*iqr[isample];

        if(is_hit && !was_hit){
	  // We just started a hit. Set the start time
	  hit.startTime=sample_time;
	  hit_charge.push_back(adc);
	  hit.SADC = adc;
	  hit.timeOverThreshold= 1;//m_downsampleFactor;
        }
        if(is_hit && was_hit){
	  hit.SADC += adc * m_downsampleFactor;
	  hit_charge.push_back(adc); 
	  hit.timeOverThreshold+=1;//m_downsampleFactor;
	}
        if(!is_hit && was_hit){
	  hit.peakCharge =  *std::max_element(hit_charge.begin(), hit_charge.end()); 
	  hit.peakTime = std::distance(hit_charge.begin(), std::max_element(hit_charge.begin(), hit_charge.end())) + hit.startTime;
	  // The hit is over. Push it to the output vector
	  hits.push_back(hit);
	  hit_charge.clear(); 
	}
        was_hit=is_hit;
    }
}


std::vector<TriggerPrimitiveFinderTool::Hit>
TriggerPrimitiveFinderPass2::findHits(const std::vector<unsigned int>& channel_numbers, 
                                      const std::vector<std::vector<short>>& collection_samples)
{
    auto hits=std::vector<TriggerPrimitiveFinderTool::Hit>();
    std::cout << "findHits called with " << collection_samples.size()
              << " channels. First chan has " << collection_samples[0].size() << " samples" << std::endl;

    for(size_t ich=0; ich<collection_samples.size(); ++ich){
        const std::vector<short>& waveformOrig=collection_samples[ich];
        std::vector<short> waveform=downSample(waveformOrig);
        std::vector<short> pedestal=findPedestal(waveform);
        std::vector<short> pedsub(waveform.size(), 0);
        for(size_t i=0; i<pedsub.size(); ++i){
            pedsub[i]=waveform[i]-pedestal[i];
        }
        std::vector<short> iqr=frugal_iqr(waveform, pedestal, m_frugalNContig);
        //skip the filtering
	//std::vector<short> filtered=filter(pedsub);
        //hitFinding(filtered, iqr, hits, channel_numbers[ich]);
	hitFinding(pedsub, iqr, hits, channel_numbers[ich]);
    }
    std::cout << "Returning " << hits.size() << " hits" << std::endl;
    std::cout << "hits/channel=" << float(hits.size())/collection_samples.size() << std::endl;
    std::cout << "hits/tick=" << float(hits.size())/collection_samples[0].size() << std::endl;
    return hits;
}

DEFINE_ART_CLASS_TOOL(TriggerPrimitiveFinderPass2)
