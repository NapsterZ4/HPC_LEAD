#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
using namespace std;

// Reads two files and compares them
int compare_files(string file_name_A, string file_name_B, double error){
	string a_line;
    string b_line;
    ifstream a_file;
    ifstream b_file;
    a_file.open(file_name_A.c_str());
    b_file.open(file_name_B.c_str());
    vector<double> a_vector;
    vector<double> b_vector;

	
    if(!a_file.is_open() || !b_file.is_open()){
        cout << "Error opening files: "<< file_name_A << " and " << file_name_B << endl;
        perror("Error open");
        // exit(EXIT_FAILURE);
	return 1; 
    }

    while(getline(a_file, a_line) && getline(b_file, b_line)) 
    {
        stringstream ss_a(a_line);
        stringstream ss_b(b_line);

        for (double i; ss_a>>i;){
            a_vector.push_back(i);
            if(ss_a.peek() == ','){
                ss_a.ignore();
            }
        }
        for (double i; ss_b>>i;){
            b_vector.push_back(i);
            if(ss_b.peek() == ','){
                ss_b.ignore();
            }
        }

        if(a_vector.size() == b_vector.size()){
            for(size_t i = 0; i<a_vector.size(); i++){
                if(fabs(a_vector[i]-b_vector[i])/fmin(fabs(a_vector[i]),fabs(b_vector[i])) > error){
                    // closing files
                    a_file.close();
                    b_file.close();
                    return 1;
                }
            }
        }else{
            // closing files
            a_file.close();
            b_file.close();
            return 1;
        }
        
    }
    a_file.close();
    b_file.close();
	return 0;
}

// Main method      
int main(int argc, char* argv[]) {
	int result;
	float error;

	// checking parameters
	if (argc != 3 && argc != 4) {
		cout<< "Parameters: <file_A> <file_B> <error>" << endl;
		return 1;
	}

	// reading error (optional)
	if(argc == 4){
		error = atof(argv[3]);
	} else {
		error = 0.2;
	}

	// comparing files, N argument is the number of lines of the files
	result = compare_files(argv[1],argv[2],error);
	cout << argv[1] << "and" << argv[2] << endl;
	if(!result)
		cout << "Tests passed" << endl;
	else
		cout << "Tests failed" << endl;
	return 0;	
}

