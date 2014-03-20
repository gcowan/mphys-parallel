#include <iostream>
void histo(){
ifstream input;
input.open("data.dat", std::ifstream::in);
TH1F * histo = new TH1F("histo","h1",50,0.0,10.0);
while(1){
    if(!input.good()) break;
    double x;
    std::string temp;
    input >> temp >> x;
    histo->Fill(x);
};
histo->Draw();
input.close();
}
