#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>



using namespace std;



string main_nuc(const string& str){
	uint Acount(0),Tcount(0);
	for(uint i(0);i<str.size();++i){
		if(str[i]=='A'){Acount++;}
		if(str[i]=='T'){Tcount++;}
	}
	if(Acount>Tcount){
		string result('A',str.size());
		return result;
	}else{
		string result('T',str.size());
		return result;
	}
	cout<<"Problem"<<endl;
	return str;
}


string clean_prefix(string str,uint min_length,uint max_missmatch,string& output){
	char c(str[0]);
	if(c!='A' and c!='T'){
		if(max_missmatch==0){
			return str;
		}else{
			return clean_prefix(str,min_length,max_missmatch-1,output);
		}
	}
	uint prefix_length(1);
	uint i(1);
	uint miss(0);
	uint con_match(0);

	for(;i< str.size();++i){
		if(str[i]==c){
			prefix_length++;
			con_match++;
		}else{
			con_match=0;
			if(++miss>max_missmatch){
				break;
			}
		}
	}
	if(con_match<5){
		prefix_length-=con_match+1;
	}
	if(prefix_length>min_length){
		output+=main_nuc(str.substr(0,i));
		str=str.substr(i);
	}
	if(str.size()==0){
		str+=output[0];
		output=output.substr(1);
	}
	return str;
}



string clean_suffix(string str,uint min_length,uint max_missmatch,string& output){
	uint j(1);
	uint suffix_length(1);
	uint miss(0);
	uint con_match(0);
	char c=(str[str.size()-1]);
	if(c!='A' and c!='T'){
		if(max_missmatch==0){
			return str;
		}else{
			return clean_suffix(str,min_length,max_missmatch-1,output);
		}
	}
	for(;j< str.size();++j){
		if(str[str.size()-1-j]==c){
			suffix_length++;
			con_match++;
		}else{
			con_match=0;
			if(++miss>max_missmatch){
				break;
			}
		}
	}
	if(con_match<5){
		suffix_length-=con_match+1;
	}
	if(suffix_length>min_length){
		output+=str.substr(str.size()-j);
		str=str.substr(0,str.size()-j);
	}

	if(str.size()==0){
		str+=output[0];
		output=output.substr(1);
	}

	return str;
}



string clean_homo(string str, uint min_length, uint max_missmatch, string& output){
	str=clean_prefix(str,min_length,max_missmatch,output);
	output+="$";
	if(str.size()<min_length){
		return str;
	}
	str=clean_suffix(str,min_length,max_missmatch,output);
	return str;
}


string getLineFasta(ifstream* in){
	string line,result;
	getline(*in,line);
	char c=in->peek();
	while(c!='>' and c!=EOF){
		getline(*in,line);
		result+=line;
		c=in->peek();
	}
	return result;
}


void clean(string& str){
	for(uint i(0); i< str.size(); ++i){
		switch(str[i]){
			case 'a':break;
			case 'A':break;
			case 'c':break;
			case 'C':break;
			case 'g':break;
			case 'G':break;
			case 't':break;
			case 'T':break;
			case 'N':break;
			default: str="";return;
		}
	}
	transform(str.begin(), str.end(), str.begin(), ::toupper);
}


int main(int argc, char ** argv){
	if(argc<2){
		cout<<"[Fasta file] (out file) (homo size) (backup file)"<<endl;
		exit(0);
	}
	string input(argv[1]);
	bool cleaning(true);
	uint min_size(0);
	string out_file("clean_reads.fa");
	string recover_file("Arecover");
	if(argc>2){
		out_file=((argv[2]));
	}
	if(argc>3){
		min_size=(stoi(argv[3]));
	}
	if(argc>4){
		recover_file=((argv[4]));
	}

	srand (time(NULL));
	string header, sequence,line;
	ifstream in(input);
	ofstream out_backup(recover_file);
	ofstream out_clean(out_file);
	string tmp_output;
	while(not in.eof()){
		getline(in,header);
		char c=in.peek();
		while(c!='>' and c!=EOF){
			getline(in,line);
			sequence+=line;
			c=in.peek();
		}
		//WE CLEAN THE SEQ
		clean(sequence);
		if(sequence.size()<5){
			sequence=clean_homo(sequence,min_size,1,tmp_output);
			out_backup<<tmp_output<<"\n";
			tmp_output="";
			out_clean<<header<<'\n'<<sequence<<"\n";
			sequence="";
		}
	}
}
