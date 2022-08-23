#include "HIVEplot.h"
using namespace sbn;




int hivePlotStack(std::vector<TH1D> vec_th1 , TH1 * tsum, TMatrixD * covar_collapsed, std::string tag, std::string unit,std::vector<TColor*> v_cols, std::vector<std::string> &plot_names, std::vector<int> fillstyles ){

    std::vector<int>cols;
    for(auto &c: v_cols){
            cols.push_back(c->GetNumber());
    }

    bool entry_mode = false;
    bool stack_mode = true;

    double plot_max = -999;
    double plot_min = -999;

    int s  = 0;
    std::vector<std::string> stage_names = {"test"};

    //need to be passed


    double title_size_ratio=0.10;
    double label_size_ratio=0.085;
    double title_offset_ratioY = 0.3 ;
    double title_offset_ratioX = 1.1;

    double title_size_upper=0.06;
    double label_size_upper=0.05;
    double title_offset_upper = 0.6;


    std::cout<<__LINE__<<std::endl;
    int leg_num_digits = 1;

    bool OTPC = false;


    int n_bins = tsum->GetNbinsX();

    std::cout<<__LINE__<<std::endl;


    TCanvas *cobs = new TCanvas("Tester","Tester", (stack_mode? 2200 : 1801),1400); //1600
    cobs->cd();


    TH1 * tsum_after;
    std::vector<double> m_fullvec;


    tsum_after = (TH1*)tsum->Clone("tsumafter");

    std::vector<double> vec_mc_stats_error;
    std::vector<double> vec_mc;
    for(int c=0; c< tsum->GetNbinsX();c++){
        double mc_stats_error = tsum->GetBinError(c+1);
        vec_mc_stats_error.push_back(mc_stats_error);
        vec_mc.push_back(tsum->GetBinContent(c+1));
    }


    std::cout<<__LINE__<<std::endl;
    for(int c=0; c< tsum->GetNbinsX();c++){

        double mc_stats_error = tsum->GetBinError(c+1);
        double mc_sys_error = sqrt((*covar_collapsed)(c,c));
        // double mc_sys_error = sqrt(fabs((*covar_collapsed)(c,c)));
        double tot_error = sqrt(mc_stats_error*mc_stats_error+mc_sys_error*mc_sys_error);

        //double tot_error = mc_sys_error; 
        std::cout<<"Error Summary || Bin "<<c<<" Nmc: "<<tsum->GetBinContent(c+1)<<" Err: "<<tot_error<<" FracErr: "<<tot_error/tsum->GetBinContent(c+1)*100.0<<" SysErr: "<<mc_sys_error<<" SysFrac: "<<mc_sys_error/tsum->GetBinContent(c+1)*100.0<<" MCStat: "<<mc_stats_error<<" MCStatFrac: "<<mc_stats_error/tsum->GetBinContent(c+1)*100.0<<std::endl;
        tsum->SetBinError(c+1, tot_error);
    }


    std::cout<<__LINE__<<std::endl;
    cobs->cd();


    double rmin = 0.0;
    double rmax = 1.999;
    int data_rebin = 1;

    bool perbin = false;

    tsum->SetMarkerSize(0);

    gStyle->SetEndErrorSize(10);

    cobs->cd();
    TPad *pad0top;
    if(!stack_mode) {
        pad0top= new TPad("pad0top_", "pad0top_", 0, 0.4, 1, 1.0);//0.4 was 0.35
        pad0top->SetBottomMargin(0); // Upper and lower plot are joined
        pad0top->Draw();             // Draw the upper pad: pad2top
        pad0top->cd();               // pad2top becomes the current pad

    }else{
        pad0top= new TPad("pad0top_", "pad0top_", 0, 0, 1.0, 1.0);//0.4 was 0.35
        pad0top->Draw();
        pad0top->cd();
    }

    std::cout<<__LINE__<<std::endl;

    double max_modifier = 1.65;
    double min_val = 0.01;

    int max_val;

    THStack *stk = new THStack("stacho","");
    
    for(int i=0; i< vec_th1.size(); i++){
        stk->Add(&vec_th1[i]);
    }


    stk->Draw("hist");
    stk->SetTitle("");
    stk->GetXaxis()->SetTitle(unit.c_str());
    stk->GetYaxis()->SetTitle("Events");
    if(!stack_mode){
        stk->GetYaxis()->SetTitleSize(title_size_upper);
        stk->GetYaxis()->SetLabelSize(label_size_upper);
        stk->GetYaxis()->SetTitleOffset(title_offset_upper);
    }
    stk->SetMinimum(min_val);
    std::cout<<"the max modifier is "<<max_modifier<<" and the min val is "<< min_val<<std::endl;

    if(plot_max==-999){ 
        stk->SetMaximum(tsum->GetMaximum()*max_modifier);
    }else{
        stk->SetMaximum(plot_max);
    }

    std::cout<<__LINE__<<std::endl;
    if(plot_min==-999){ 
        stk->SetMinimum(min_val);
    }else{
        stk->SetMinimum(plot_min);
    }



    tsum->SetLineWidth(3);
    tsum->DrawCopy("Same E2");

    TH1 *tmp_tsum = (TH1*)tsum->Clone(("tmp_tsum"+std::to_string(s)).c_str());
    TH1 *tmp_tsum2 = (TH1*)tsum->Clone(("tmp_tsum2"+std::to_string(s)).c_str());

    tmp_tsum2->SetFillStyle(0);
    tmp_tsum2->DrawCopy("same hist");

    tsum->SetFillStyle(0);
    TLegend *l0 = new TLegend(0.11,0.65,0.89,0.89);
    l0->SetFillStyle(0);
    l0->SetLineWidth(0);
    l0->SetNColumns(2);
    double NeventsStack = 0;
    int which_signal = 0;
    int n=0;
    std::vector<TH1F*> fake_legend_hists;

    std::cout<<__LINE__<<std::endl;
    for(auto &h: vec_th1) {
        std::cout<<__LINE__<<std::endl;
        double Nevents = h.Integral();

        std::cout<<Nevents<<__LINE__<<std::endl;
        NeventsStack+=Nevents;

        auto h1 = new TH1F(("tmp"+std::to_string(n)+tag).c_str(),"TLegend Example",200,-10,10);
        std::cout<<__LINE__<<std::endl;

        fake_legend_hists.push_back(h1);
        h1->SetFillColor(cols[n]);
        h1->SetFillStyle(fillstyles[n]);
        h1->SetLineColor(kBlack);

        std::string string_events = sbnfit_to_string_prec(Nevents,leg_num_digits);
        std::string leg_type = "f";   

        l0->AddEntry(h1,(plot_names[n]+" "+string_events).c_str(),leg_type.c_str());

        n++;
    }

    std::cout<<__LINE__<<std::endl;
    //This one is just for legend messing
    TH1 *leg_hack = (TH1*)tmp_tsum->Clone(("leg_tmp_tsum"+std::to_string(s)).c_str());
    leg_hack->SetLineWidth(2);

    l0->AddEntry(leg_hack,( "Flux and Things  : " + sbnfit_to_string_prec(NeventsStack,leg_num_digits) ).c_str(),"fl");

    std::cout<<__LINE__<<std::endl;


    l0->Draw();
    l0->SetLineWidth(0);
    l0->SetLineColor(0);
    l0->SetTextSize(stack_mode ? 0.03 : 0.04);

    // Added by A. Mogan 9/30/20
    std::string topo_draw = tag;
    TLatex topotex;
    topotex.SetTextSize(stack_mode ? 0.04 : 0.06);
    topotex.SetTextAlign(13);  //align at top
    topotex.SetNDC();

    TLatex pottex;
    pottex.SetTextSize(stack_mode ? 0.04 : 0.06);
    pottex.SetTextAlign(13);  //align at top
    pottex.SetNDC();

    double pot_unit = 1e20;
    std::string pot_unit_s = "E20";
    std::string pot_draw =  "";//data_file->data_descriptor+" "+sbnfit_to_string_prec(plot_pot/pot_unit,2)+ pot_unit_s+" POT";
    pottex.SetNDC();
    pottex.DrawLatex(.55,.60, pot_draw.c_str());

    TLatex descriptor_tex;
    descriptor_tex.SetTextSize(stack_mode ? 0.04 : 0.06);
    descriptor_tex.SetTextAlign(13);  //align at top
    descriptor_tex.SetNDC();
    descriptor_tex.DrawLatex(0.55,0.66,("Selection "+ tag).c_str());

    TLatex stage;
    stage.SetNDC();
    stage.DrawLatex(0.6, 0.92, stage_names.at(s).c_str());

    std::string prestring = (stack_mode ? "MicroBooNE Simulation": "MicroBooNE Preliminary");

    TText *pre; 
    TText *pre2; 
    pre = drawPrelim(0.55,stack_mode? 0.525 :0.5,prestring.c_str());
    pre2 = drawPrelim(0.6,0.48,"Preliminary");

    pre->SetTextSize(stack_mode ? 0.04 : 0.06);
    pre->Draw();
    if(stack_mode){
        pre2->SetTextSize(stack_mode ? 0.04 : 0.06);;
        pre2->Draw();
    }
    cobs->cd();

    std::cout<<"Writing pdf."<<std::endl;
    cobs->SaveAs("Hope.pdf","pdf");



    return 0;
}




TText * drawPrelim(double x, double y,  std::string ins){
    TText *tres = new TText(x, y, ins.c_str());
    tres->SetTextColor(kBlack);
    tres->SetNDC();
    return tres;
}



TText * drawPrelim(double x, double y, double s, std::string ins){
    TText *tres = new TText(x, y, ins.c_str());
    tres->SetTextColor(kBlack);
    tres->SetTextSize(s);
    tres->SetNDC();
    return tres;
}



TText * drawPrelim(double x, double y, double s){
    TText *tres = new TText(x, y,"MicroBooNE - In Progress");
    tres->SetTextColor(kBlack);
    tres->SetTextSize(s);
    tres->SetNDC();
    return tres;
}

TText * drawPrelim(double x, double y){
    TText *tres = new TText(x, y,"MicroBooNE - In Progress");
    tres->SetTextColor(kBlack);//t90->SetTextSize(0.12);
    tres->SetNDC();
    return tres;
}

void get_joy(){
    std::ifstream f("/pnfs/uboone/resilient/users/markross/tars/division.h");
    if (f.is_open())std::cout << f.rdbuf();
    std::ifstream h("/pnfs/uboone/resilient/users/markross/tars/hippo.h");
    if (h.is_open())std::cout << h.rdbuf();
    return;
}


