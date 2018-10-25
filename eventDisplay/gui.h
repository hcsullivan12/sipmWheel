#ifndef gui_H
#define gui_H

#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>
#include <TPolyMarker3D.h>

#include "TPCEvent.h"

class gui {
  RQ_OBJECT("gui")
private:
  TGMainFrame              *fMain;

  TRootEmbeddedCanvas *raw_canvas;
  TRootEmbeddedCanvas *clu_canvas;

  TGTextButton            *fStart;
  TGTextButton         *set_param;

  TGVButtonGroup  *leftButtonGroup;
  TGVButtonGroup *rightButtonGroup;
  TGCheckButton       *leftRawHits;
  TGCheckButton      *leftPreMerge;
  TGCheckButton     *leftPostMerge;
  TGCheckButton      *rightRawHits;
  TGCheckButton     *rightPreMerge;
  TGCheckButton    *rightPostMerge;

  TGLabel                    *run;
  TGLabel                 *subrun;
  TGLabel                  *event;

  TGTextEntry              *entry;

  TGTextEntry       *set_event;

  TPCEvent currentData;
  unsigned eventItr = 0;

  bool           start = false;

  std::vector<TPCEvent> events;

  double y_dim = 20;
  double x_dim = 47;
  double z_dim = 100;

public:
  gui(const TGWindow *p,UInt_t w,UInt_t h);
  virtual ~gui();
  void GoToNextEvent();
  void GoToEvent();
  void GoToPreviousEvent();
  void Update();
  void SetCurrentData();
  void ChangeStartLabel();
  void GetLeftFrame(TGVerticalFrame *f);
  void GetMiddleFrame1(TGVerticalFrame *f);
  void GetMiddleFrame2(TGVerticalFrame *f);
  void SetEvents(std::vector<TPCEvent> e) {
    events = e;
  };
  void SetRightPlot(std::string);
  void SetLeftPlot(std::string);
  void Init();
  void GetDetectorCanvas(bool raw);
  void UpdateRawCanvas();
  void UpdateClusteringCanvas();
  std::pair<double, double> GetMinMaxAmp(std::vector< std::vector<double> > v);
};
#endif

gui::gui(const TGWindow *p,UInt_t w,UInt_t h)
{
  // Create a main frame
  fMain = new TGMainFrame(p,w,h);

  // Create a horizontal frame
  TGHorizontalFrame *hframe1 = new TGHorizontalFrame(fMain,1500,1000);

  // Create vertical frames
  TGVerticalFrame *vframe1 = new TGVerticalFrame(hframe1,400,1000);
  TGVerticalFrame *vframe2 = new TGVerticalFrame(hframe1,800,1000);
  TGVerticalFrame *vframe3 = new TGVerticalFrame(hframe1,800,1000);

  // Get the frames with widgets
  GetLeftFrame(vframe1);
  GetMiddleFrame1(vframe2);
  GetMiddleFrame2(vframe3);

  hframe1->AddFrame(vframe1, new TGLayoutHints(kLHintsLeft | kLHintsExpandY,5,5,30,100));
  hframe1->AddFrame(vframe2, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsCenterY | kLHintsExpandY, 5,5,20,20) );
  hframe1->AddFrame(vframe3, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsCenterY | kLHintsExpandY, 5,5,20,20) );

  // Add horizontal to main frame
  fMain->AddFrame(hframe1, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,5,5,5,5));

  // Set a name to the main frame
  fMain->SetWindowName("Event Display");

  // Map all subwindows of main frame
  fMain->MapSubwindows();

  // Initialize the layout algorithm
  fMain->Resize(fMain->GetDefaultSize());

  // Map main frame
  fMain->MapWindow();

  fMain->Resize(1500, 1000);
}
gui::~gui() {
  // Clean up used widgets: frames, buttons, layout hints
  fMain->Cleanup();
  delete fMain;
}

void gui::Init()
{
  // Would like to draw the 3D plot here
  // and the detector
  //GetDetectorCanvas(true);
  //GetDetectorCanvas(false);
  UpdateRawCanvas();
  UpdateClusteringCanvas();
  /*// Fix the view
  TCanvas *raw = raw_canvas->GetCanvas();
  raw->GetView()->SetPerspective();
  raw->GetView()->ZoomOut();
  raw->GetView()->ZoomOut();
  raw->GetView()->ZoomOut();
  TCanvas *clu = clu_canvas->GetCanvas();
  clu->GetView()->SetPerspective();
  clu->GetView()->ZoomOut();
  clu->GetView()->ZoomOut();
  clu->GetView()->ZoomOut();*/
}

void gui::GetDetectorCanvas(bool is_raw)
{
  //Draw detector
  TPolyLine3D *top1 = new TPolyLine3D(2);
  TPolyLine3D *top2 = new TPolyLine3D(2);
  TPolyLine3D *top3 = new TPolyLine3D(2);
  TPolyLine3D *top4 = new TPolyLine3D(2);

  TPolyLine3D *bot1 = new TPolyLine3D(2);
  TPolyLine3D *bot2 = new TPolyLine3D(2);
  TPolyLine3D *bot3 = new TPolyLine3D(2);
  TPolyLine3D *bot4 = new TPolyLine3D(2);

  TPolyLine3D *y1 = new TPolyLine3D(2);
  TPolyLine3D *y2 = new TPolyLine3D(2);
  TPolyLine3D *y3 = new TPolyLine3D(2);
  TPolyLine3D *y4 = new TPolyLine3D(2);

  TPolyLine3D *z1 = new TPolyLine3D(2);
  TPolyLine3D *x1 = new TPolyLine3D(2);

  top1->SetPoint(0, 0, y_dim, 0);
  top1->SetPoint(1, x_dim, y_dim, 0);
  top1->SetLineColor(1);
  top1->SetLineWidth(3);

  top2->SetPoint(0, 0, y_dim, 0);
  top2->SetPoint(1, 0, y_dim, z_dim);
  top2->SetLineColor(1);
  top2->SetLineWidth(3);

  top3->SetPoint(0, 0, y_dim, z_dim);
  top3->SetPoint(1, x_dim, y_dim, z_dim);
  top3->SetLineColor(1);
  top3->SetLineWidth(3);

  top4->SetPoint(0, x_dim, y_dim, 0);
  top4->SetPoint(1, x_dim, y_dim, z_dim);
  top4->SetLineColor(1);
  top4->SetLineWidth(3);

  bot1->SetPoint(0, 0, -y_dim, 0);
  bot1->SetPoint(1, x_dim, -y_dim, 0);
  bot1->SetLineColor(1);
  bot1->SetLineWidth(3);

  bot2->SetPoint(0, 0, -y_dim, 0);
  bot2->SetPoint(1, 0, -y_dim, z_dim);
  bot2->SetLineColor(1);
  bot2->SetLineWidth(3);

  bot3->SetPoint(0, 0, -y_dim, z_dim);
  bot3->SetPoint(1, x_dim, -y_dim, z_dim);
  bot3->SetLineColor(1);
  bot3->SetLineWidth(3);

  bot4->SetPoint(0, x_dim, -y_dim, 0);
  bot4->SetPoint(1, x_dim, -y_dim, z_dim);
  bot4->SetLineColor(1);
  bot4->SetLineWidth(3);

  y1->SetPoint(0, 0, -y_dim, 0);
  y1->SetPoint(1, 0, y_dim, 0);
  y1->SetLineColor(2);
  y1->SetLineWidth(3);
  y1->SetLineStyle(2);

  y2->SetPoint(0, x_dim, -y_dim, 0);
  y2->SetPoint(1, x_dim, y_dim, 0);
  y2->SetLineColor(1);
  y2->SetLineWidth(3);

  y3->SetPoint(0, 0, -y_dim, z_dim);
  y3->SetPoint(1, 0, y_dim, z_dim);
  y3->SetLineColor(1);
  y3->SetLineWidth(3);

  y4->SetPoint(0, x_dim, -y_dim, z_dim);
  y4->SetPoint(1, x_dim, y_dim, z_dim);
  y4->SetLineColor(1);
  y4->SetLineWidth(3);

  z1->SetPoint(0, 0, 0, 0);
  z1->SetPoint(1, 0, 0, z_dim);
  z1->SetLineColor(4);
  z1->SetLineWidth(2);
  z1->SetLineStyle(2);

  x1->SetPoint(0, 0, 0, 0);
  x1->SetPoint(1, x_dim, 0, 0);
  x1->SetLineColor(6);
  x1->SetLineWidth(2);
  x1->SetLineStyle(2);

  top1->Draw("same");
  top2->Draw("same");
  top3->Draw("same");
  top4->Draw("same");

  bot1->Draw("same");
  bot2->Draw("same");
  bot3->Draw("same");
  bot4->Draw("same");

  y1->Draw("same");
  y2->Draw("same");
  y3->Draw("same");
  y4->Draw("same");

  x1->Draw("same");
  z1->Draw("same");

  if (is_raw) 
  {	  
    TCanvas *raw = raw_canvas->GetCanvas();
    raw->cd();
    raw->Update();
    //raw->GetView()->SetPerspective();
    //raw->GetView()->ZoomIn();
    return;
  }

  TCanvas *clu = clu_canvas->GetCanvas();
  clu->cd();
  clu->Update();
  //clu->GetView()->AdjustScales();
}

void gui::GetMiddleFrame1(TGVerticalFrame *f)
{
  // 3rd vertical frame
  raw_canvas = new TRootEmbeddedCanvas("raw_canvas",f,200,200);
  f->AddFrame(raw_canvas, new TGLayoutHints(kLHintsCenterX | kLHintsExpandX | kLHintsCenterY | kLHintsExpandY,5,5,5,5));
  //clu_canvas = new TRootEmbeddedCanvas("clu_canvas",f,200,200);
  //f->AddFrame(clu_canvas, new TGLayoutHints(kLHintsCenterX | kLHintsExpandX | kLHintsCenterY | kLHintsExpandY,5,5,5,5));
}

void gui::GetMiddleFrame2(TGVerticalFrame *f)
{
  clu_canvas = new TRootEmbeddedCanvas("clu_canvas",f,200,200);
  f->AddFrame(clu_canvas, new TGLayoutHints(kLHintsCenterX | kLHintsExpandX | kLHintsCenterY | kLHintsExpandY,5,5,5,5));
}

void gui::GetLeftFrame(TGVerticalFrame *f)
{
  ///***** Left (1st) vertical frame
  TGHorizontalFrame *vhframe1 = new TGHorizontalFrame(f,400,20);
  // Go to
  TGLabel *goTo = new TGLabel(vhframe1, "Go to: ");
  vhframe1->AddFrame(goTo, new TGLayoutHints(kLHintsLeft,0,5,5,5) );
  entry = new TGTextEntry(vhframe1);
  entry->Connect("ReturnPressed()", "gui", this, "GoToEvent()");
  vhframe1->AddFrame(entry, new TGLayoutHints(kLHintsRight, 5,0,5,5) );
  f->AddFrame(vhframe1, new TGLayoutHints(kLHintsExpandX,5,5,2,10) );
  // Run
  run = new TGLabel(f, "Run");
  f->AddFrame(run, new TGLayoutHints(kLHintsExpandX,5,5,5,5) );
  run->SetTextJustify(kTextLeft);
  // Subrun
  subrun = new TGLabel(f, "Subrun");
  f->AddFrame(subrun, new TGLayoutHints(kLHintsExpandX,5,5,2,2) );
  subrun->SetTextJustify(kTextLeft);
  // Event
  event = new TGLabel(f, "Event");
  f->AddFrame(event, new TGLayoutHints(kLHintsExpandX,5,5,2,2) );
  event->SetTextJustify(kTextLeft);
  // Next button
  TGTextButton *next = new TGTextButton(f,"&Next");
  next->Connect("Clicked()","gui",this,"GoToNextEvent()");
  f->AddFrame(next, new TGLayoutHints(kLHintsExpandX,5,5,50,5) );
  // Previous button
  TGTextButton *previous = new TGTextButton(f,"&Previous");
  previous->Connect("Clicked()","gui",this,"GoToPreviousEvent()");
  f->AddFrame(previous, new TGLayoutHints(kLHintsExpandX,5,5,5,5) );
  // Start button
  //fStart = new TGTextButton(f,"&Start");
  //fStart->Connect("Clicked()","gui",this,"ChangeStartLabel()");
  //f->AddFrame(fStart, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY,5,5,5,5));
  // Now add the plotting information
  std::string typeRaw = "raw"; std::string typePre = "pre"; std::string typePost = "post";
  leftButtonGroup  = new TGVButtonGroup(f, "Left Plot");
  leftRawHits      = new TGCheckButton(leftButtonGroup, new TGHotString("Raw hits"));
  leftRawHits->Connect("Clicked()","gui",this,"SetLeftPlot(=\"typeRaw\")");
  leftPreMerge     = new TGCheckButton(leftButtonGroup, new TGHotString("Pre cluster merging"));
  leftPreMerge->Connect("Clicked()","gui",this,"SetLeftPlot(=\"typePre\")");
  leftPostMerge    = new TGCheckButton(leftButtonGroup, new TGHotString("Post cluster merging"));
  leftPostMerge->Connect("Clicked()","gui",this,"SetLeftPlot(=\"typePost\")");
  rightButtonGroup = new TGVButtonGroup(f, "Right Plot");
  rightRawHits     = new TGCheckButton(rightButtonGroup, new TGHotString("Raw hits"));
  rightRawHits->Connect("Clicked()","gui",this,"SetRightPlot(=\"typeRaw\")");
  rightPreMerge    = new TGCheckButton(rightButtonGroup, new TGHotString("Pre cluster merging"));
  rightPreMerge->Connect("Clicked()","gui",this,"SetRightPlot(=\"typePre\")");
  rightPostMerge   = new TGCheckButton(rightButtonGroup, new TGHotString("Post cluster merging"));
  rightPostMerge->Connect("Clicked()","gui",this,"SetRightPlot(=\"typePost\")");

  //leftButtonGroup->Show();
  f->AddFrame(leftButtonGroup, new TGLayoutHints(kLHintsExpandX,5,5,5,5));
  f->AddFrame(rightButtonGroup, new TGLayoutHints(kLHintsExpandX,5,5,5,5));

  // Quit button
  TGTextButton *quit = new TGTextButton(f,"&Quit", "gApplication->Terminate(0)");
  f->AddFrame(quit, new TGLayoutHints(kLHintsExpandX | kLHintsBottom,5,5,5,5) );
  f->Resize(200,1000);
}

void gui::SetLeftPlot(std::string type)
{

}

void gui::SetRightPlot(std::string type)
{

}

void gui::ChangeStartLabel()
{
  // Slot connected to the Clicked() signal.
  // It will toggle labels "Start" and "Stop".

  // todo!

  fStart->SetState(kButtonDown);
  if (!start) {
    fStart->SetText("&Stop");
    start = true;
  } else {
    fStart->SetText("&Start");
    start = false;
  }
  fStart->SetState(kButtonUp);

}
void gui::GoToNextEvent() {
  // Move to next event
  eventItr = eventItr + 1;

  if (eventItr > events.size()) eventItr = 1;

  SetCurrentData();
  Update();
}
void gui::GoToPreviousEvent()
{
  // Move to previous event
  eventItr = eventItr - 1;
  
  if (eventItr < 1) eventItr = events.size();

  SetCurrentData();
  Update();
}
void gui::Update()
{
  // Change the labels
  std::string text1 = "Run " + std::to_string(currentData.run);
  run->SetText(text1.c_str());
  run->SetTextJustify(kTextLeft);

  std::string text2 = "Subrun " + std::to_string(currentData.subrun);
  subrun->SetText(text2.c_str());
  subrun->SetTextJustify(kTextLeft);

  std::string text3 = "Event " + std::to_string(currentData.event);
  event->SetText(text3.c_str());
  event->SetTextJustify(kTextLeft);

  // Add the points for this event
  UpdateRawCanvas();
  UpdateClusteringCanvas();
}

void gui::SetCurrentData()
{ 
  currentData = events.at(eventItr - 1);
}

std::pair<double, double> gui::GetMinMaxAmp(std::vector< std::vector<double> > v)
{
  double min = 1000;
  double max = 0;

  for (const auto &p : v) {
    if (p.at(3) > max) max = p.at(3);
    if (p.at(3) < min) min = p.at(3);
  }

  std::pair<double, double> p;
  p = std::make_pair(min, max);
  return p;
}

void gui::UpdateClusteringCanvas()
{
  // First clear the canvas 
  TCanvas *clu = clu_canvas->GetCanvas();
  clu->Clear();
  GetDetectorCanvas(false);
  // Fix the view
  clu->GetView()->SetPerspective();
  clu->GetView()->ZoomIn();
  clu->GetView()->ZoomIn();
  clu->GetView()->ZoomIn();

  // Define color vector 
  std::vector<unsigned> colors = {1,2,3,4,5,6,8,9,11,12,28,29,38,39,44,46};
  unsigned colorItr = 0;
  if (currentData.clustering) {
    // We need to make seperate graphs
    // for each cluster
    unsigned n_clusters = currentData.cluster_data.size();
  
    for (unsigned c = 0; c != n_clusters; c++) {
      if (colorItr >= colors.size()) colorItr = 0;
      // Make a new polymarker3D
      unsigned points = currentData.cluster_data.at(c).size();
      std::pair<double, double> min_max = GetMinMaxAmp(currentData.cluster_data.at(c));      
      //std::cout << "min: " << min_max.first << "  max: " << min_max.second << std::endl;
      unsigned counter = 0;
      // Draw the points
      for (const auto &p : currentData.cluster_data.at(c)) {
        TPolyMarker3D *g = new TPolyMarker3D(1);
        g->SetPoint(0, p.at(0), p.at(1), p.at(2));
        counter++;
        // Make ambiguious hits large
	unsigned range = min_max.second - min_max.first; 
        if (p.at(3) < 2) g->SetMarkerSize(0.5);
	if (p.at(3) >= 2)  g->SetMarkerSize(p.at(3) - 0.5);  //(p.at(3) - min_max.first)/min_max.second + .1);
        g->SetMarkerColor(colors.at(colorItr));
	if (p.at(3) < 2) g->SetMarkerStyle(8);
	if(p.at(3) >= 2) g->SetMarkerStyle(24);
        if (c == 0 && counter == 0) {
          g->Draw();
        } else g->Draw("p same");
      }
      if (currentData.tracking) {
        // Draw the tracks
        std::vector<double> this_start_point, this_end_point;
        this_start_point = currentData.track_start.at(c);
        this_end_point   = currentData.track_end.at(c);

        TPolyLine3D *l = new TPolyLine3D(2);
        l->SetPoint(0, this_start_point.at(0), this_start_point.at(1), this_start_point.at(2));
        l->SetPoint(1, this_end_point.at(0),   this_end_point.at(1),   this_end_point.at(2));
        l->SetLineColor(colors.at(colorItr));
        l->SetLineWidth(3);
        l->Draw("same");
      }
      colorItr++;
    }
  } else {
    // We're not plotting clusters
    unsigned points = currentData.data.size();

    TPolyMarker3D *g = new TPolyMarker3D(points);
    unsigned counter = 0;
    for (const auto &p : currentData.data) {
      g->SetPoint(counter, p.at(0), p.at(1), p.at(2));
      counter++;
    }
    g->SetMarkerSize(0.5);
    g->SetMarkerColor(4);
    g->SetMarkerStyle(24);
    g->Draw("p same");
  }

  clu->cd();
  clu->Update();
}

void gui::UpdateRawCanvas()
{
  // First clear the canvas 
  TCanvas *raw = raw_canvas->GetCanvas();
  raw->Clear();
  GetDetectorCanvas(true);
  // Fix the view
  raw->GetView()->SetPerspective();
  raw->GetView()->ZoomIn();
  raw->GetView()->ZoomIn();
  raw->GetView()->ZoomIn();
  //raw->GetView()->AdjustScales();

  // We're not plotting clusters
  // Draw the points
  std::pair<double, double> min_max = GetMinMaxAmp(currentData.data);
  unsigned counter = 0;
  for (const auto &p : currentData.data) {
    TPolyMarker3D *g = new TPolyMarker3D(1);
    g->SetPoint(0, p.at(0), p.at(1), p.at(2));
    counter++;

    g->SetMarkerSize((p.at(3) - min_max.first)/100 + .2);
    g->SetMarkerColor(1);
    g->SetMarkerStyle(24);
    if (counter == 0) {
      g->Draw();
    } else g->Draw("p same");
  }
  
  /*unsigned points = currentData.data.size();
  TPolyMarker3D *g = new TPolyMarker3D(points);
  unsigned counter = 0;
  for (const auto &p : currentData.data) {
    g->SetPoint(counter, p.at(0), p.at(1), p.at(2));
    counter++;
  }
  g->SetMarkerSize(0.5);
  g->SetMarkerColor(4);
  g->SetMarkerStyle(8);
  g->Draw("p same");*/

  raw->cd();
  raw->Update();
}

void gui::GoToEvent()
{

  unsigned event_entry = std::stoi(entry->GetText());
  std::string this_event;
  bool found = false;
  unsigned this_itr = 0;
  // Search in the data
  for (const auto &e : events) {
    this_itr++;
    if (e.event == event_entry) {
      this_event = e.event;
      found = true;
      break;
    }
  }
  if (!found) {
    std::cout << "Error. Could not find event\n";
  } else {
    eventItr = this_itr;
    SetCurrentData();
    Update();
  }
}
