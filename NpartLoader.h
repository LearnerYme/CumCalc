#include <fstream>
#include <string>
#include <vector>
#include <iostream>

class NpartLoader{

    private:
        static const int nCent = 9;
        int edge[nCent];

    public:
        NpartLoader(){
            std::ifstream fin;
            fin.open("Npart.txt");
            std::string str;
            int i = 0;
            while(std::getline(fin, str)){
                edge[i] = std::atoi(str.c_str());
                i += 1;
                if (i > nCent){
                    std::cout << "[Warning] In CentDefinition class, the target centrality edge is out of range " << i << " of " << nCent << ".\n";
                }
            }
            fin.close();
        }
        ~NpartLoader(){}

        int* GetArray() { return edge; }
};
