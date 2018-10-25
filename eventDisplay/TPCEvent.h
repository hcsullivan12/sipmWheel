#ifndef TPCEVENT_HH
#define TPCEVENT_HH

struct TPCEvent {
  
  public:
    float run;
    float subrun;
    float event;
    bool clustering = false;
    bool tracking = false;
    std::vector< std::vector<double> > data;
    std::vector< std::vector< std::vector<double> > > cluster_data; 
    std::vector< std::vector<double> > track_start;
    std::vector< std::vector<double> > track_end;
};

#endif
