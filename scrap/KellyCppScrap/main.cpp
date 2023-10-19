#include <iostream>
#include <vector>
#include <matrix.h>
#include <KellyCppScrap.h>

using namespace std;

int main()
{
    unsigned int NSpec = 2;
    unsigned int NComp = 2;
    vector<float> CConc[NComp] = {
        1.00748E-07,// H(+)
        3.90357E-08 // CO3(2-)
    };
    vector<float> K[NSpec] = {
      2.133045e+10, //K_HCO3 = 10^10.329
      4.852885e+16  //K_H2CO3 = 10^16.686
    }
    signed int StoicM[NSpec][NComp] = {
        {1, 1}, // HCO3 = 1*H + 1*CO3
        {2, 1}  // H2CO3 = 2*H + 1*CO3
    };
    Matrix<signed int> Stoic(NSpec, NComp, StoicM[0], NSpec * NComp);
    vector<float> SConc;

    try {
        SConc = CppCalcSpecConc(Stoic, CConc);
        for (iSpec = 0; iSpec < NSpec; iSpec++) {
            std::cout <<  SConc[iSpec] << "  ";
        }
    } catch (MatrixException& e) {
        std::cerr << e.message() << std::endl;
        return e.errorCode();
    }

    return 0;
}
