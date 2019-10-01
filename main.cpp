#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <math.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <utility>
#include <bits/stdc++.h>
#include <sys/time.h>

using namespace std; 

struct timeval start_time, end_time;

class block
{
public:
	int name;
	double Height;
	double Width;
	double x_coord;                  // Lower X-coordinate of the Block
	double y_coord;	              // Lower Y-coordinate of the Block
	int pos_position;             // Position of the Block in the Positive Sequence
	int neg_position;             // Position of the Block in the Negative Sequence
	vector<block*> fanIn_HCG;     // Fan-in for the Block in the Horizontal Graph
	vector<block*> fanOut_HCG;    // Fan-Out for the Block in the Horizontal Graph
	vector<block*> fanIn_VCG;     // Fan-in for the Block in the Vertical Graph
	vector<block*> fanOut_VCG;    // Fan-Out for the Block in the Vertical Graph
	int ip_cnt_hcg;
	int ip_cnt_vcg;
public: void Block_initialize(int Block_Name, double w, double h)
 {
	name = Block_Name;
	Height = h;	
	Width = w;	
 }
};

///////////////////////// Vectors and Vector Lists //////////////////////////////////////
vector<block> B;
vector<block*> B_ptr;
vector<list<block*> > vector_blocks;
vector<block*> pos_seq;
vector<block*> neg_seq;

/////////////////////// Queues used for HCG and VCGs ////////////////////////////////////
queue<block*> queue_HCG;
queue<block*> queue_VCG;


/////////////////////////// Functions Declaration ///////////////////////////////////////
void swap_positive(int x, int y);
void swap_negative(int x, int y);
void swap_both(int x, int y);
double calculate_wirelength();
void calculate_HCG_VCG();
void calculate_fanin_out();
double calculate_Area();
bool AcceptedMove(double Delta,double Temperature);
///////////////////////////////////////////////////////////////////////////////////////// 

//////////////////////// Global Variable Declaration ////////////////////////////////////
int Width_Initial;
int Height_Initial;
int Width_Final;
int Height_Final;
/////////////////////////////////////////////////////////////////////////////////////////

class Net
{
public:
	string name;
	int degree;
};


int main(int argc, char **argv)
{
////////////////////// Time Calculation ////////////////////////////

    double start_count, end_count;
    double elapsed_time;        
    gettimeofday(&start_time,NULL);

///////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////
//clock_t tStart = clock();
//std::srand(std::time(0));
//////////////////////////////////////////////////////////////////
srand(time(NULL));
if(argc!=3)
 {
 cout<<"Provide valid Input Arguments"<<endl;
 cout<<"Example: ./parser <filename> -a" << endl;
 return 1;
 }

string filename=argc[argv-2];
string Annealing_Type=argc[argv-1];

int block_size=0, net_size=0, out=0;
string Block_Name,n_name;
int name;
double w,h;
size_t st,end,nt,del,sz;


ifstream input_file(filename.c_str()); 
string line_n;
if(!input_file.is_open())
{
	cout << "File could not be opened.\n";
	return 1;
}
while(getline(input_file, line_n))
{

	st = line_n.find_first_not_of(" ",st);
	string block_size_str = line_n.substr(0,st);

	block_size = atoi(block_size_str.c_str());

//cout<<"Block Size: " << block_size << endl;

/////////////////////////////////// Parsing the Bench File /////////////////////////////////////


////////////////////// First Parse the Blocks and store them in B vector of block //////////////////////

	while(getline(input_file, line_n))
	{
		block b;
		if((nt=line_n.find("Nets")) == string::npos)
		{	
			
			st = line_n.find_first_not_of(" ", 0);
			end = line_n.find_first_of(" ", st);
			Block_Name = line_n.substr(st,end);
			name = atoi(Block_Name.c_str());
			st = line_n.find_first_not_of(" ", end);
			end = line_n.find_first_of(" ", st);
			string s_w = line_n.substr(st,end-st);
			w = atoi(s_w.c_str());
			st = line_n.find_first_not_of(" ", end);
			end = line_n.find_last_of(" ", st);
			string s_h = line_n.substr(st,end-st);
			h = atoi(s_h.c_str());
			//cout << "name " << name << " Width " << w << " Height " << h << endl;
			b.Block_initialize(name,w,h);		
			
		}
		else
		{ 
			getline(input_file,line_n);	
			st = line_n.find_first_not_of(" ",st);
			string net_size_str = line_n.substr(0,st);
			net_size = atoi(net_size_str.c_str());
			//cout<<" net_size: " << net_size << endl;
			break; 
		}
	B.push_back(b);	
	}

////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////// Next moving on to Parsing the Nets //////////////////////////////////////////////	

	while(getline(input_file, line_n))
	{	
		list<block*> local_list;
		vector<block>::iterator it;
		
		st = line_n.find_first_not_of(" ", 0);
		end = line_n.find_first_of(" ", st);
		n_name = line_n.substr(st,end);
		name = atoi(n_name.c_str());
		do
		{
			st = line_n.find_first_not_of(" ", end);
			end = line_n.find_first_of(" ", st);
			string B_ptr = line_n.substr(st,end-st);
			out = atoi(B_ptr.c_str());
			for(it = B.begin(); it!= B.end(); it++)
			{
				if(it->name == out)
					break;
			}
			local_list.push_back(&(*it));
		//cout << "Out: " << out << endl;
		}while(end != string::npos);

		vector_blocks.push_back(local_list);                             ////////// Each line pushed to the vector_blocks as a pointer
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////// Storing in Pointers !! /////////////////////////////////////

for(int i=0;i<B.size();i++)
	B_ptr.push_back(&B[i]);


/////////////////////////// Check for an Initial Sequence Hardcoded by the User ////////////////////////////////
/*
B_ptr[0]->pos_position = 0;
B_ptr[0]->neg_position = 7;
B_ptr[1]->pos_position = 4;
B_ptr[1]->neg_position = 3;
B_ptr[2]->pos_position = 6;
B_ptr[2]->neg_position = 5;
B_ptr[3]->pos_position = 2;
B_ptr[3]->neg_position = 1;
B_ptr[4]->pos_position = 3;
B_ptr[4]->neg_position = 4;
B_ptr[5]->pos_position = 5;
B_ptr[5]->neg_position = 6;
B_ptr[6]->pos_position = 1;
B_ptr[6]->neg_position = 2;
B_ptr[7]->pos_position = 7;
B_ptr[7]->neg_position = 0;
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////// Print out the Blocks from list of Pointer Blocks //////////////////////////////////

list<block*>::iterator it;
for(int i=0; i<vector_blocks.size(); i++)
{
	for(it = vector_blocks[i].begin(); it != vector_blocks[i].end(); it++)
	{
			block* temp = *it;
			//cout << "Block name: " << temp->name << " Width: " << temp->Width << " Height: " << temp->Height << endl; 
	}
	//cout << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////// First Define a Random Positive and Negative Sequence ///////////////////////////////
////////////////////////// The Sequences are stored in pos_seq and neg_seq vectors ////////////////////////////

for(int i=0; i<B_ptr.size(); i++)
{
	B_ptr[i]->pos_position = i;
	B_ptr[i]->neg_position = i;
	pos_seq.push_back(B_ptr[i]);
	neg_seq.push_back(B_ptr[i]);	
}

///////////////////////////  Initial Run //////////////////////////////////
calculate_fanin_out();
calculate_HCG_VCG();
//cout<< "Chip Width: " << Width_Final << " Chip Height: " << Height_Final << endl;
//cout << endl <<endl;

double Initial_Wire_Length = calculate_wirelength();
//cout << "Initial Wire Length: " << Initial_Wire_Length << endl;

////////////////////////////////////////////////////////////////////////////
Width_Initial = Width_Final;
Height_Initial = Height_Final;
////////////////////////////////////////////////////////////////////////////

double Initial_Area = calculate_Area();
//cout << "Initial Area of Chip: " << Initial_Area << endl;

////////////////////////////////////////////////////////////////////////////



//////////////////////// Simulated Annealing ///////////////////////////////

double current_wirelength = Initial_Wire_Length;
double next_wirelength = 0;
double Delta_wirelength = 0;

double current_Area = Initial_Area;
double next_Area = 0;
double Delta_Area = 0;

double Delta_Total = 0;

double initial_temperature = 40000;
double freezing_temperature = 1;
double Temperature = initial_temperature;

int random1,random2,move;
double Temperature_Step = 0.95;

bool obj = 0;

/////////////////////////////////////////////////////////////////////////////
string var7 = argc[argv-2];
string constant_var7 = "_Results.csv";
string var8;
if(Annealing_Type == "-a")
	var8 = "a";
if(Annealing_Type == "-w")
	var8 = "w";
if(Annealing_Type == "-c")
	var8 = "c"; 
////////////////////////////////////////////////////////////////////////////
string res1_filename = var7+var8+constant_var7;

ofstream res7(res1_filename.c_str());
res7<<"Temperature"<<"\t"<<"Accepted Moves"<<"\t"<<"Rejected Moves\t"<<"HPWL\t"<<"Area"<<endl;
ofstream res1("Results.csv");
res1<<"Temperature"<<"\t"<<"Accepted Moves"<<"\t"<<"Rejected Moves\t"<<"HPWL\t"<<"Area"<<endl;

while(Temperature>freezing_temperature)
{
	int accepted = 0;
	int rejected = 0;
	for(int i=0;i<100;i++)
	{
		random1 = rand()%B_ptr.size();
		random2 = rand()%B_ptr.size();
		move = rand()%3;

		if(move==0)
			swap_positive(random1,random2);
		else if(move == 1)
			swap_negative(random1,random2);
		else
			swap_both(random1,random2);
		
		/////////// For WireLength ////////////////////
		if(Annealing_Type == "-w")
		{
		//cout << "Doing for Wirelength" << endl;
		next_wirelength = calculate_wirelength();
		Delta_wirelength = next_wirelength - current_wirelength;
		//current_Area = calculate_Area();
		bool a = AcceptedMove(Delta_wirelength,Temperature);

		if(a)
		{
			current_wirelength = next_wirelength;
			accepted++;
		}
		else
		{
			rejected++;
			if(move==0)
				swap_positive(random1,random2);
			else if(move == 1)
				swap_negative(random1,random2);
			else
				swap_both(random1,random2);
		}
		current_Area = calculate_Area();
		}
		/////////////////////////////////////////////

		///////////// For Area //////////////////////
		else if(Annealing_Type == "-a")
		{
		//cout << "Doing for Area" << endl;
		next_Area = calculate_Area();
		Delta_Area = next_Area - current_Area;
		//current_wirelength = calculate_wirelength();
		bool b = AcceptedMove(Delta_Area,Temperature);

		if(b)
		{
			current_Area = next_Area;
			accepted++;
		}
		else
		{
			rejected++;
			if(move==0)
				swap_positive(random1,random2);
			else if(move == 1)
				swap_negative(random1,random2);
			else
				swap_both(random1,random2);
		}
		current_wirelength = calculate_wirelength();
		}
		////////////////////////////////////////////

		//////////// For Both Area and Wirelength ///////////////////////////
		else
		{
		//cout << "Doing for Both" << endl;
		next_wirelength = calculate_wirelength();
		Delta_wirelength = next_wirelength - current_wirelength;

		next_Area = calculate_Area();
		Delta_Area = next_Area - current_Area;

		double alpha = 0.2;
		Delta_Total = alpha*(Delta_Area) + (1-alpha)*(Delta_wirelength);
		
		bool c = AcceptedMove(Delta_Total,Temperature);

		if(c)
		{
			current_Area = next_Area;
			current_wirelength = next_wirelength;
			accepted++;
		}
		else
		{
			rejected++;
			if(move==0)
				swap_positive(random1,random2);
			else if(move == 1)
				swap_negative(random1,random2);
			else
				swap_both(random1,random2);
		}
		}
		/////////////////////////////////////////////////////////////////////
	}
//cout << Temperature << endl;


//place<<"Final Placement of Gates"<<endl;
res1<<Temperature<<"\t"<<accepted<<"\t"<<rejected << "\t"<< current_wirelength<< "\t"<< current_Area <<endl;
res7<<Temperature<<"\t"<<accepted<<"\t"<<rejected << "\t"<< current_wirelength<< "\t"<< current_Area <<endl;
Temperature = Temperature*Temperature_Step;
}

//cout << "Initial Wire Length: " << Initial_Wire_Length << endl;
//cout << "Final Wirelength: " << current_wirelength << endl;


//cout << "Initial Area: " << Initial_Area << endl;
//cout << "Final Area: " << current_Area << endl;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////// ### Write Everything to a File //////////////////////////////////////////////
string var1 = argc[argv-2];
string constant_val = "_Sheeparamatti_Arshath.out2";
string var2; 
//cout << Annealing_Type << endl;
if(Annealing_Type == "-a")
	var2 = "a";
if(Annealing_Type == "-w")
	var2 = "w";
if(Annealing_Type == "-c")
	var2 = "c";
string outfile_name=var1+constant_val+var2;

ofstream res2(outfile_name.c_str());
res2 << "Results: " << endl;
res2 << "Initial Wirelength of the Placement is: " << Initial_Wire_Length << endl;
res2 << "Initial Width of the Placement is: " << Width_Initial << endl;
res2 << "Initial Height of the Placement is: " << Height_Initial << endl;
res2 << "Initial Area of the Placement is: " << Initial_Area << endl;
res2 << endl;
res2 << "------------------------------------------------------------------------------" << endl;
res2 << "The Total HPWL of the Placement is: " << current_wirelength << endl;
res2 << "The width of the Placement is: " << Width_Final << endl;
res2 << "The height of the Placement is: " << Height_Final << endl;
res2 << "The Total Area of Placement is: " << current_Area << endl;
res2 << endl;
res2 << "------------------------------------------------------------------------------" << endl;
res2 << "Final Block Placements with x and y coordinates: " << endl;

for(int i=0;i<B_ptr.size(); i++)
{
	res2<< "Block Name: " << (*B_ptr[i]).name << "\t" <<  "Block x-coodinate: " << (*B_ptr[i]).x_coord << "\t" <<  "Block y-coordinate: " << (*B_ptr[i]).y_coord << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
for(int i=0; i<B_ptr.size(); i++)
{
	for(int j=0; j< (*B_ptr[i]).fanIn_VCG.size(); j++)
		cout<< "\nIn VCG Block name: " << (*B_ptr[i]).name << " Block inputs name: " << (*B_ptr[i]).fanIn_VCG[j]->name<< endl;
	for(int j=0; j<B_ptr[i]->fanOut_VCG.size(); j++)
		cout<< "In VCG Block name: " << (*B_ptr[i]).name << " Block outputs name: " << (*B_ptr[i]).fanOut_VCG[j]->name<< endl;
	for(int j=0; j< (*B_ptr[i]).fanIn_HCG.size(); j++)
		cout<< "In HCG Block name: " << (*B_ptr[i]).name << " Block inputs name: " << (*B_ptr[i]).fanIn_HCG[j]->name<< endl;
	for(int j=0; j<B_ptr[i]->fanOut_HCG.size(); j++)
		cout<< "In HCG Block name: " << (*B_ptr[i]).name << " Block outputs name: " << (*B_ptr[i]).fanOut_HCG[j]->name<< endl;
}
*/

//cout << endl << endl;
//for(int i=0; i<B_ptr.size(); i++)
//{
//	cout << B_ptr[i]->name << "\t" << B_ptr[i]->x_coord << "\t" << B_ptr[i]->Width << "\t"<<endl;
//}


//swap_positive(0,2);
//swap_both(3,5);
//Final_Wire_Length = calculate_wirelength();
//cout << "Final Wire Length: " << Final_Wire_Length << endl;

///////////////////////////////////// Print them out and check ////////////////////////////////////////////////
/*
cout << endl << endl;
for(int i=0; i<B_ptr.size(); i++)
{
	cout << B_ptr[i]->name << "\t" << B_ptr[i]->x_coord << "\t" << B_ptr[i]->Width << "\t"<<endl;
}
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//printf("Time Taken: %.3fs\n",(double)(clock()-tStart)/CLOCKS_PER_SEC);


    gettimeofday(&end_time,NULL);
    start_count = (double) start_time.tv_sec + 1.e-6 * (double) start_time.tv_usec;
    end_count = (double) end_time.tv_sec + 1.e-6 * (double) end_time.tv_usec;
    elapsed_time = (end_count - start_count);
    printf("The total elapsed time is:  %f seconds\n",elapsed_time);
    res2 << "\nEXECUTION TIME IS : " << "\t" << elapsed_time << endl;


return 0;

}



////////////////////////////////////// Functions Definitions ////////////////////////////////////////////////////

////////////////////////// ### Function to calculate Fan-in and Fan-out /////////////////////////////////////////
//////////////////////////       for each individual Block in the FP   //////////////////////////////////////////

void calculate_fanin_out()
{
	vector<block*>::iterator iter;
	iter=B_ptr.begin();
	vector<block*>::iterator iter1;
	iter1=B_ptr.begin();

	for(iter=B_ptr.begin(); iter != B_ptr.end(); iter++)
	{
		(*iter)->fanIn_HCG.clear();
		(*iter)->fanOut_HCG.clear();
		(*iter)->fanIn_VCG.clear();
		(*iter)->fanOut_VCG.clear();
	}


	for(iter=B_ptr.begin();iter!=B_ptr.end();iter++)
	{
		for(iter1=B_ptr.begin();iter1!=B_ptr.end();iter1++)
		{
	   		if((*iter)->name != (*iter1)->name)
	    		{
				if((*iter)->pos_position > (*iter1)->pos_position && (*iter)->neg_position < (*iter1)->neg_position)
				{
					(*iter)->fanOut_VCG.push_back(&(*(*iter1)));
					(*iter1)->fanIn_VCG.push_back(&(*(*iter)));
					//cout << (*iter)->name << endl;
                        		//cout<<"inside"<<endl;
				}

				if((*iter)->pos_position < (*iter1)->pos_position && (*iter)->neg_position < (*iter1)->neg_position)
				{
					(*iter)->fanOut_HCG.push_back(&(*(*iter1)));
					(*iter1)->fanIn_HCG.push_back(&(*(*iter)));
					//cout << (*iter)->name << endl;
					//cout<<"Entering Once"<<endl;
				}
	    		}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////// ### Function to calculate HCG and the VCG //////////////////////////////////////
////////////////////////////////        for each iteration of Annealing   ///////////////////////////////////////

void calculate_HCG_VCG()
{
	Width_Final = 0;
	Height_Final = 0;
	
	///////////////////////////////////// Traverse the HCG and VCG Graphs //////////////////////////////////////////
	for(int k=0; k<B_ptr.size(); k++)
	{
		B_ptr[k]->ip_cnt_hcg = (*B_ptr[k]).fanIn_HCG.size();
		B_ptr[k]->ip_cnt_vcg = (*B_ptr[k]).fanIn_VCG.size();
		B_ptr[k]->x_coord = 0;
		B_ptr[k]->y_coord = 0;
		if((*B_ptr[k]).fanIn_HCG.size() == 0)
		{
			(*B_ptr[k]).x_coord = 0;
			//cout << "Block: " << (*B_ptr[k]).name << "\tX-coordinate: " << (*B_ptr[k]).x_coord << endl;
			queue_HCG.push(B_ptr[k]);
		}

		if((*B_ptr[k]).fanIn_VCG.size() == 0)
		{
			(*B_ptr[k]).y_coord = 0;
			//cout << "Block: " << (*B_ptr[k]).name << "\tY-coordinate: " << (*B_ptr[k]).y_coord << endl;
			queue_VCG.push(B_ptr[k]);
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////// Traverse the VCG and assign the y-positions ////////////////////////////////////
	//cout << queue_VCG.size() << endl; 
	while(!queue_VCG.empty())
	{
		//cout << "VCG Queue" << endl;
		block* b = queue_VCG.front();
		for(int i=0; i < (*b).fanOut_VCG.size(); i++)
		{
			//cout << "Entering For" << endl;
			block* Out = (*b).fanOut_VCG[i];
			Out->ip_cnt_vcg--;
			////////////////////////// Check for the Maximum Distance ///////////////////////////
			if(b->y_coord+b->Height > Out->y_coord)
				Out->y_coord = b->y_coord+b->Height;

			// Could be commented back in later
			//if(Height_Final < Out->y_coord + Out->Height)
				//Height_Final = Out->y_coord + Out->Height;
				//cout << "Entering If" <<endl;
			/////////////////////////////////////////////////////////////////////////////////////		
			if(Out->ip_cnt_vcg == 0)
			{
				//Out->y_coord = b->y_coord+b->Height;
				//cout << "Entering here" << endl;
				queue_VCG.push(Out);
				//cout << "Block Name: " << Out->name << "\tBlock y-coordinate: " << Out->y_coord << endl;
			}

		}
		queue_VCG.pop();
	}


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////// Traverse the HCG and assign the x-positions ///////////////////////////////////

	while(!queue_HCG.empty())
	{
		//cout << "HCG Queue" << endl;
		block* b = queue_HCG.front();
		for(int i=0; i < (*b).fanOut_HCG.size(); i++)
		{
			block* Out = (*b).fanOut_HCG[i];
			Out->ip_cnt_hcg--;
			////////////////////////// Check for the Maximum Distance ///////////////////////////
			if(b->x_coord + b->Width > Out->x_coord)
				Out->x_coord = b->x_coord + b->Width;
				//cout << "Doin it! X-Coordinate" << endl;

			/////////////////////////////////////////////////////////////////////////////////////
			// Could be commented back in later			
			//if(Width_Final < Out->x_coord + Out->Width)
			//	Width_Final = Out->x_coord + Out->Width;
			/////////////////////////////////////////////////////////////////////////////////////

			/////////////////////////////////////////////////////////////////////////////////////		
			if(Out->ip_cnt_hcg == 0)
			{
				//Out->x_coord = b->x_coord+b->Width;
				queue_HCG.push(Out);
				//cout << "Block Name: " << Out->name << "\tBlock x-coordinate: " << Out->x_coord << endl;
			}

		}
		queue_HCG.pop();
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////// Final Width and Height Calculations for the Chip ///////////////////////////////
	
	for(int k=0;k<B_ptr.size();k++)
	{
		if(Width_Final < (*B_ptr[k]).x_coord + (*B_ptr[k]).Width)
			Width_Final = (*B_ptr[k]).x_coord + (*B_ptr[k]).Width;
		if(Height_Final < (*B_ptr[k]).y_coord + (*B_ptr[k]).Height)
			Height_Final = (*B_ptr[k]).y_coord + (*B_ptr[k]).Height;
	}	

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////


}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////// ### Function to Swap 2 elements of the Positive Sequence /////////////////////////////

void swap_positive(int x, int y)
{
	//cout << "Blocks being swapped are: " << B_ptr[x]->name <<" " << B_ptr[y]->name << endl;
	swap(B_ptr[x]->pos_position,B_ptr[y]->pos_position);
	//cout << endl << endl;
	
	//cout<< "Whatever follows is after the Swap" << endl;
	calculate_fanin_out();
	calculate_HCG_VCG();
	//cout<< "Chip Width: " << Width_Final << " Chip Height: " << Height_Final << endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////// ### Function to Swap 2 elements of the Negative Sequence /////////////////////////////

void swap_negative(int x, int y)
{
	//cout << "Blocks being swapped are: " << B_ptr[x]->name <<" " << B_ptr[y]->name << endl;
	swap(B_ptr[x]->neg_position,B_ptr[y]->neg_position);
	
	
	calculate_fanin_out();
	calculate_HCG_VCG();
	//cout<< "Chip Width: " << Width_Final << " Chip Height: " << Height_Final << endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////// ### Function to Swap 2 elements of Both Sequences ////////////////////////////////////

void swap_both(int x, int y)
{
	//cout << "Blocks being swapped are: " << B_ptr[x]->name <<" " << B_ptr[y]->name << endl;
	swap(B_ptr[x]->pos_position,B_ptr[y]->pos_position);
	swap(B_ptr[x]->neg_position,B_ptr[y]->neg_position);


	calculate_fanin_out();
	calculate_HCG_VCG();
	//cout<< "Chip Width: " << Width_Final << " Chip Height: " << Height_Final << endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////// ### Function to Calculate the Wire Length //////////////////////////////////////
double calculate_wirelength()
{
	double hpwl = 0;
	vector<double> temp1;
	vector<double> temp2;
	for(int i=0; i<vector_blocks.size(); i++)
	{
		list<block*> u = vector_blocks[i];
		list<block*>::iterator iter;
		//vector<double> temp1;
		//vector<double> temp2;
		vector<double>::iterator maximum_x,maximum_y,minimum_x,minimum_y;

		for(iter = u.begin(); iter!=u.end();iter++)
		{
			block *b = *iter;
			temp1.push_back(b->x_coord+(b->Width/2));
			temp2.push_back(b->y_coord+(b->Height/2));
		}
		
		maximum_x = max_element(temp1.begin(),temp1.end());
		maximum_y = max_element(temp2.begin(),temp2.end());
		minimum_y = min_element(temp2.begin(),temp2.end());
		minimum_x = min_element(temp1.begin(),temp1.end());

		double inter = (*maximum_x-*minimum_x)+(*maximum_y-*minimum_y);
		hpwl = hpwl + inter;
		temp1.clear();
		temp2.clear();
	}
	//temp1.clear();
	//temp2.clear();
	return hpwl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////// ### Function to Calculate the Intermediate Areas //////////////////////////////////
double calculate_Area()
{
	double area=0;
	return Width_Final*Height_Final;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////// ### Accept Move /////////////////////////////////////
bool AcceptedMove(double Delta,double Temperature)
{
	if(Delta < 0) return 1;
	double boltzman_const = exp(((-1)*(Delta))/((Temperature*0.95)));
	double random_compare =((double)rand()/(double)RAND_MAX);
	if(random_compare<boltzman_const)
		return 1;
	else return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////

