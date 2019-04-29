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
		string result(str.size(),'A');
		return result;
	}else{
		string result(str.size(),'T');
		return result;
	}
	cout<<"Problem"<<endl;
	return str;
}



pair<string,string> clean_prefix(string str,uint min_length,uint max_missmatch,string output){
	//~ return  {str,output};
	//~ cout<<"cp"<<endl;
//~ cout<<min_length<<endl;
	char c(str[0]);
	if(c!='A' and c!='T'){
		if(max_missmatch==0){
			return {str,output};
		}else{
			auto rec(clean_prefix(str.substr(1),min_length,max_missmatch-1,output));
			if(rec.first.size()==str.size()-1){
				return {str,output};
			}else{
				return {rec.first,c+rec.second};
			}
		}
	}
	uint prefix_length(1);
	uint i(1);
	uint miss(0);
	uint con_match(1);

	for(;i< str.size();++i){
		prefix_length++;
		if(str[i]==c){
			con_match++;
		}else{
			if(++miss>max_missmatch){
				break;
			}
			con_match=1;
		}
	}
	if(con_match<5){
		//~ cout<<"remove first macthes"<<endl;
		prefix_length-=(con_match+1);
	}else{
		prefix_length--;
		//~ cout<<"keep first macthes"<<endl;
	}
	if(prefix_length>min_length){
		output+=main_nuc(str.substr(0,prefix_length));
		str=str.substr(prefix_length);
	}
	if(str.size()==0){
		str+=output[0];
		output=output.substr(1);
	}
	return {str,output};
}


uint count_upper_case(const string& str){
	int res(0);
	for(uint i(0); i< str.size(); ++i){
		switch(str[i]){
			//~ case 'a':break;
			case 'A':++res;break;
			//~ case 'c':break;
			case 'C':++res;break;
			//~ case 'g':break;
			case 'G':++res;break;
			//~ case 't':break;
			case 'T':++res;break;
			//~ default: return res;
		}
	}
	return res;
}



pair<string,string> clean_suffix(string str, uint min_length, uint max_missmatch,  string output){
	//~ return  {str,output};
	//~ cout<<"cs"<<endl;
	uint j(1);
	uint suffix_length(1);
	uint miss(0);
	uint con_match(1);
	char c=(str[str.size()-1]);
	if(c!='A' and c!='T'){
		if(max_missmatch==0){
			return {str,output};
		}else{
			auto rec(clean_suffix(str.substr(0,str.size()-1),min_length,max_missmatch-1,output));
			if(rec.first.size()==str.size()-1){
				return {str,output};
			}else{
				return {rec.first,rec.second+c};
			}
		}
	}
	for(;j< str.size();++j){
		suffix_length++;
		if(str[str.size()-1-j]==c){
			con_match++;
		}else{
			if(++miss>max_missmatch){
				break;
			}
			con_match=1;
		}
	}
	//~ cout<<"suffix length"<<endl;
	//~ cout<<suffix_length<<endl;
	if(con_match<5){
		//~ cout<<"remove last macthes"<<endl;
		suffix_length-=(con_match+1);
	}else{
		suffix_length--;
		//~ cout<<"keep last macthes"<<endl;
	}
	//~ cout<<"suffix length"<<endl;
	//~ cout<<suffix_length<<endl;

	if(suffix_length>min_length){
		output+=main_nuc(str.substr(str.size()-suffix_length));
		str=str.substr(0,str.size()-suffix_length);
	}
	//~ cout<<str<<endl;

	if(str.size()==0){
		str+=output[0];
		output=output.substr(1);
	}
	//~ cout<<"output"<<output<<endl;
	return  {str,output};
}


pair<string,string> clean_prefix2(const string& str, uint min_length, uint max_missmatch){
	if(str.size()<min_length){
		return {str,""};
	}
	uint ca(0),cc(0),cg(0),ct(0);
	for(uint i(0);i<min_length;++i){
		switch(str[i]){
			case 'A':++ca;break;
			case 'C':++cc;break;
			case 'G':++cg;break;
			default:++ct;break;
		}
		if(cc+cg>max_missmatch){
			return {str,""};
		}
	}
	if(ca>ct){
		if (ca<min_length-max_missmatch){
			return {str,""};
		}
	}else{
		if (ct<min_length-max_missmatch){
			return {str,""};
		}
	}
	uint nuc_to_remove(0);
	for(uint i(0);i+min_length<str.size();++i){
		switch(str[i]){
			case 'A':--ca;break;
			case 'C':--cc;break;
			case 'G':--cg;break;
			default:--ct;break;
		}
		switch(str[i+min_length]){
			case 'A':++ca;break;
			case 'C':++cc;break;
			case 'G':++cg;break;
			default:++ct;break;
		}
		if(ca>ct){
			if (ca<min_length-max_missmatch){
				break;
			}
		}else{
			if (ct<min_length-max_missmatch){
				break;
			}
		}
		nuc_to_remove++;
	}
	return{str.substr(nuc_to_remove+min_length),str.substr(0,nuc_to_remove)};
}



pair<string,string> clean_homo(string& str, uint min_length, uint max_missmatch){
	if(str.size()<min_length){
		return {str, ""};
	}
	string output;
	auto pair=clean_prefix(str,min_length,max_missmatch,output);
	pair.second+="$";
	if(pair.first.size()<min_length){
		return pair;
	}
	auto pair2=clean_suffix(pair.first,min_length,max_missmatch,pair.second);
	return pair2;
}



pair<string,string> clean_homo2(string& str, uint min_length, uint max_missmatch){
	auto pair=clean_prefix2(str,min_length,max_missmatch);
	reverse(pair.first.begin(),pair.first.end());
	auto pair2=clean_prefix2(pair.first,min_length,max_missmatch);
	reverse(pair2.first.begin(),pair2.first.end());
	return {pair2.first,main_nuc(pair.second)+"$"+main_nuc(pair2.second)};
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
	uint min_size(10);
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
	uint i(0);
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
		if(sequence.size()>5){
			auto pair=clean_homo2(sequence,min_size,2);
			out_backup<<pair.second<<"\n";
			tmp_output="";
			out_clean<<header<<'\n'<<pair.first<<"\n";
		}
		sequence="";
	}
}
