#include <fstream>
#include <string>
#include <iostream>

class CentDefinition{

    private:
        static const int nCent = 9;
        int edge[nCent];
        int is_set;

    public:
        CentDefinition(){
            for (int i=0; i<nCent; i++){
                edge[i] = 0;
            }
            is_set = 0;
        }
        ~CentDefinition(){}

        void read_edge(const char* path){
            std::ifstream fin;
            fin.open(path);
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
            is_set = 1;
            std::cout << "[LOG] Set centrality edge done.\n";
            std::cout << "[LOG] Edge: ";
            for (int i=0; i<nCent; i++){
                std::cout << edge[i];
                if (i != nCent-1){
                    std::cout << ", ";
                } else {
                    std::cout << ".\n";
                }
            }
        }

        int get_cent(int mult){
            if (!is_set){
                std::cout << "[Warning] You should read centrality bin edge first.\n";
                return -1;
            }
            for (int i=0; i<nCent; i++){
                if (mult > edge[i]){
                    return i;
                }
            }
            return -1;
        }
};
