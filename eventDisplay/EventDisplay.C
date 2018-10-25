#include "gui.h"
#include "FileReader.h"
#include "TPCEvent.h"

void EventDisplay(std::string input_file) 
{ 
  // First need to read in the data
  FileReader fr(input_file);
  fr.ReadFile();

  std::vector<TPCEvent> events = fr.GetEvents();
  //fr.PrintClusterData();

  // Open the GUI...
  auto g = new gui(gClient->GetRoot(),1500,1000);
  g->SetEvents(events);
  g->Init();
}

