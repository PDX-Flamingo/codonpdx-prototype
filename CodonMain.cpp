#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <string>

using namespace std;

class List {
public:
    Node* first;
    Node* last;

    List()
    {
        first = NULL;
        last = NULL;
        listLength = 0;
    };


    long listLength;
    void AddLinkToBack(hostOrganism* ptr) {
        Node* newNode = new Node;
        newNode->data = ptr;
        Node * current = first;

        if (listLength == 0) {
            first = newNode;
            last = newNode;
            newNode->next = NULL;
            newNode->previous = NULL;
        }
        else {
            while (current->next !=NULL) {
                current = current->next;
            }

            current->next = newNode;
            newNode->previous = current;
            newNode->next = NULL;
            last = newNode;

        }
        listLength++;
    }

};


struct Node {
    Node* next;
    Node* previous;
    hostOrganism* data;

    Node() {
        next = NULL;
        previous = NULL;
        data = NULL;
    }
};
// global constants
const int MAX_CHAR = 101;
const int AA = 64;
const int MAX_HOSTS = 100;
const int codonLength = 2000;
List* llptr;





struct hostOrganism {
    char hostNumber[MAX_CHAR];
    char hostName[MAX_CHAR];
    int codons;
    double codonUsage[AA];
    double usagePercentage;
};





struct inputGenome {
    double codonUsage[AA];
    double proteinTotals[AA];
};

//prototypes
void openFile(ifstream& inFile, ofstream& outFile, int& place);
void convertDNAtoRNA(char genome[], int sizeOfArray);
void sortResults(int sizeOfArray, hostOrganism addHost[]);
void sortGenome(char genome[], int sizeOfArray, inputGenome newGenome[], int inputProtein[]);
void incrementInputCodon(inputGenome newGenome[], char firstCodon, char secondCodon, char thirdCodon, int inputProtein[]);
void calculatePercentUsage(inputGenome newGenome[], int inputProtein[]);


int main()
{
    char stop;
    int place = 0;
    ifstream inFile;
    ofstream outFile;
    llptr = new List();

    //hostOrganism addHost[MAX_HOSTS];
    char genome[codonLength] = {'\0'};
    //char userChoice = ' ';
    int inputProtein[21] = {0}; // R,L,S,T,P,A,G,V,K,N,Q,H,E,D,Y,C,F,I,M,W,STOP

    openFile(inFile, outFile, place);

//while(tolower(userChoice) != x) {
    inputGenome newGenome[1];
//cout << "Press 'x' to exit, or any other key to continue: ");
//cin >> userChoice;
    cout << "\nEnter the sequence to be analyzed: " << endl;
    cin.getline(genome, MAX_CHAR, '\n');
    int sizeOfArray = strlen(genome);
//cout << sizeOfArray << endl;

    convertDNAtoRNA(genome, sizeOfArray);
    sortGenome(genome, sizeOfArray, newGenome, inputProtein);
    calculatePercentUsage(newGenome, inputProtein);

//}

    cin >> stop;
    return 0;
}


// **********************

void openFile(ifstream& inFile, ofstream& outFile,hostOrganism addHost[], int& place) {

    char stop;
    int tempCodonNumber;
    double tempCodonUsage[AA];

// opens the file test.txt

    inFile.open("test.txt");
    if(!inFile.is_open()) // if file is unable to open
    {
        cout << "Unable to save new file." << endl;
    }

    // reads from the existing file
    while (!inFile.eof())	{
        string temp;
        hostOrganism* newData = new hostOrganism;

        inFile.getline(newData->hostNumber, MAX_CHAR, ':');
        inFile.getline(newData->hostName, MAX_CHAR, ':');
        inFile >> tempCodonNumber;
        newData->codons = tempCodonNumber;

        for(int i = 0; i < AA; i++) {
            inFile >> newData->codonUsage[i];
            //tempCodonUsage[i];
        }

        llptr->AddLinkToBack(newData);

        /*
        		cout << addHost[place].hostNumber << " " << addHost[place].hostName << " " << addHost[place].codons << endl;
        		for(int i = 0; i < AA; i++) {
        			cout << tempCodonUsage[i] << "\n";
        			cout << addHost[place].codonUsage[i] << " ";
        		} */
    }
    //cout << addHost[10].hostNumber << " " << addHost[10].hostName << " " << addHost[10].codons << endl;
    //cin >> stop;
    inFile.close();
    return;
}


// ******************************

void convertDNAtoRNA(char genome[], int sizeOfArray) {

    for(int i = 0; i < sizeOfArray; i++) {
        if(genome[i] == 'T') {
            genome [i] = 'U';
        }
    }

    return;
}


//*******************

void sortResults(int sizeOfArray, hostOrganism addHost[]) {

    hostOrganism temp;

    for(int i = 0; i < sizeOfArray; i++) {
        for(int j = i+1; j < sizeOfArray; j++) {
            if(addHost[i].usagePercentage > addHost[j].usagePercentage) {
                temp = addHost[i];
                addHost[i] = addHost[j];
                addHost[j] = temp;

            }
        }
    }
    return;
}


//***********************

void sortGenome(char genome[], int sizeOfArray, inputGenome newGenome[], int inputProtein[]) {

    char firstCodon = ' ';
    char secondCodon = ' ';
    char thirdCodon = ' ';

    for(int i = 0; i < sizeOfArray; i + 3) {
        for(int j = 1; j < (i + 3); j + 3) {
            for(int k = 2; k < (i + 3); k + 3) {

                genome[i] = firstCodon;
                genome[j] = secondCodon;
                genome[k] = thirdCodon;

                incrementInputCodon(newGenome, firstCodon, secondCodon, thirdCodon, inputProtein);

            }
        }
    }


    return;
}


//************************

void incrementInputCodon(inputGenome newGenome[], char firstCodon, char secondCodon, char thirdCodon, int inputProtein[]) {

    int R = 0;
    int L = 0;
    int	S = 0;
    int T = 0;
    int P = 0;
    int A = 0;
    int G = 0;
    int V = 0;
    int K = 0;
    int N = 0;
    int Q = 0;
    int H = 0;
    int E = 0;
    int D = 0;
    int Y = 0;
    int C = 0;
    int F = 0;
    int I = 0;
    int M = 0;
    int W = 0;
    int STOP = 0;


    /*
    CGA = 0,R  CGC = 1,R  CGG = 2,R  CGU = 3,R  AGA = 4,R  AGG = 5,R  CUA = 6,L
    CUC = 7,L  CUG = 8,L  CUU = 9,L  UUA = 10,L UUG = 11,L UCA = 12,S UCC = 13,S
    UCG = 14,S UCU = 15,S AGC = 16,S AGU = 17,S ACA = 18,T ACC = 19,T ACG = 20,T
    ACU = 21,T CCA = 22,P CCC = 23,P CCG = 24,P CCU = 25,P GCA = 26,A GCC = 27,A
    GCG = 28,A GCU = 29,A GGA = 30,G GGC = 31,G GGG = 32,G GGU = 33,G GUA = 34,V
    GUC = 35,V GUG = 36,V GUU = 37,V AAA = 38,K AAG = 39,K AAC = 40,N AAU = 41,N
    CAA = 42,Q CAG = 43,Q CAC = 44,H CAU = 45,H GAA = 46,E GAG = 47,E GAC = 48,D
    GAU = 49,D UAC = 50,Y UAU = 51,Y UGC = 52,C UGU = 53,C UUC = 54,F UUU = 55,F
    AUA = 56,I AUC = 57,I AUU = 58,I AUG = 59,M UGG = 60,W UAA = 61,STOP
    UAG =62,STOP UGA = 63,STOP

    */

    if(firstCodon == 'C' && secondCodon == 'G' && thirdCodon == 'A') {
        newGenome[0].codonUsage[0]= newGenome[0].codonUsage[0] + 1;
        R++;
    }

    if(firstCodon == 'C' && secondCodon == 'G' && thirdCodon == 'C') {
        newGenome[0].codonUsage[1] = newGenome[0].codonUsage[1] + 1;
        R++;
    }

    if(firstCodon == 'C' && secondCodon == 'G' && thirdCodon == 'G') {
        newGenome[0].codonUsage[2] = newGenome[0].codonUsage[2] + 1;
        R++;
    }

    if(firstCodon == 'C' && secondCodon == 'G' && thirdCodon == 'U') {
        newGenome[0].codonUsage[3]++;
        R++;
    }

    if(firstCodon == 'A' && secondCodon == 'G' && thirdCodon == 'A') {
        newGenome[0].codonUsage[4]++;
        R++;
    }

    if(firstCodon == 'A' && secondCodon == 'G' && thirdCodon == 'G') {
        newGenome[0].codonUsage[5]++;
        R++;
    }

    if(firstCodon == 'C' && secondCodon == 'U' && thirdCodon == 'A') {
        newGenome[0].codonUsage[6]++;
        L++;
    }

    if(firstCodon == 'C' && secondCodon == 'U' && thirdCodon == 'C') {
        newGenome[0].codonUsage[7]++;
        L++;
    }

    if(firstCodon == 'C' && secondCodon == 'U' && thirdCodon == 'G') {
        newGenome[0].codonUsage[8]++;
        L++;
    }

    if(firstCodon == 'C' && secondCodon == 'U' && thirdCodon == 'U') {
        newGenome[0].codonUsage[9]++;
        L++;
    }

    if(firstCodon == 'U' && secondCodon == 'U' && thirdCodon == 'A') {
        newGenome[0].codonUsage[10]++;
        L++;
    }

    if(firstCodon == 'U' && secondCodon == 'U' && thirdCodon == 'G') {
        newGenome[0].codonUsage[11]++;
        L++;
    }

    if(firstCodon == 'U' && secondCodon == 'C' && thirdCodon == 'A') {
        newGenome[0].codonUsage[12]++;
        S++;
    }

    if(firstCodon == 'U' && secondCodon == 'C' && thirdCodon == 'C') {
        newGenome[0].codonUsage[13]++;
        S++;
    }

    if(firstCodon == 'U' && secondCodon == 'C' && thirdCodon == 'G') {
        newGenome[0].codonUsage[14]++;
        S++;
    }

    if(firstCodon == 'U' && secondCodon == 'C' && thirdCodon == 'U') {
        newGenome[0].codonUsage[15]++;
        S++;
    }

    if(firstCodon == 'A' && secondCodon == 'G' && thirdCodon == 'C') {
        newGenome[0].codonUsage[16]++;
        S++;
    }

    if(firstCodon == 'A' && secondCodon == 'G' && thirdCodon == 'U') {
        newGenome[0].codonUsage[17]++;
        S++;
    }

    if(firstCodon == 'A' && secondCodon == 'C' && thirdCodon == 'A') {
        newGenome[0].codonUsage[18]++;
        T++;
    }

    if(firstCodon == 'A' && secondCodon == 'C' && thirdCodon == 'C') {
        newGenome[0].codonUsage[19]++;
        T++;
    }

    if(firstCodon == 'A' && secondCodon == 'C' && thirdCodon == 'G') {
        newGenome[0].codonUsage[20]++;
        T++;
    }

    if(firstCodon == 'A' && secondCodon == 'C' && thirdCodon == 'U') {
        newGenome[0].codonUsage[21]++;
        T++;
    }

    if(firstCodon == 'C' && secondCodon == 'C' && thirdCodon == 'A') {
        newGenome[0].codonUsage[22]++;
        P++;
    }

    if(firstCodon == 'C' && secondCodon == 'C' && thirdCodon == 'C') {
        newGenome[0].codonUsage[23]++;
        P++;
    }

    if(firstCodon == 'C' && secondCodon == 'C' && thirdCodon == 'G') {
        newGenome[0].codonUsage[24]++;
        P++;
    }

    if(firstCodon == 'C' && secondCodon == 'C' && thirdCodon == 'U') {
        newGenome[0].codonUsage[25]++;
        P++;
    }

    if(firstCodon == 'G' && secondCodon == 'C' && thirdCodon == 'A') {
        newGenome[0].codonUsage[26]++;
        A++;
    }

    if(firstCodon == 'G' && secondCodon == 'C' && thirdCodon == 'C') {
        newGenome[0].codonUsage[27]++;
        A++;
    }

    if(firstCodon == 'G' && secondCodon == 'C' && thirdCodon == 'G') {
        newGenome[0].codonUsage[28]++;
        A++;
    }

    if(firstCodon == 'G' && secondCodon == 'C' && thirdCodon == 'U') {
        newGenome[0].codonUsage[29]++;
        A++;
    }

    if(firstCodon == 'G' && secondCodon == 'G' && thirdCodon == 'A') {
        newGenome[0].codonUsage[30]++;
        G++;
    }

    if(firstCodon == 'G' && secondCodon == 'G' && thirdCodon == 'C') {
        newGenome[0].codonUsage[31]++;
        G++;
    }

    if(firstCodon == 'G' && secondCodon == 'G' && thirdCodon == 'G') {
        newGenome[0].codonUsage[32]++;
        G++;
    }

    if(firstCodon == 'G' && secondCodon == 'G' && thirdCodon == 'U') {
        newGenome[0].codonUsage[33]++;
        G++;
    }

    if(firstCodon == 'G' && secondCodon == 'U' && thirdCodon == 'A') {
        newGenome[0].codonUsage[34]++;
        V++;
    }

    if(firstCodon == 'G' && secondCodon == 'U' && thirdCodon == 'C') {
        newGenome[0].codonUsage[35]++;
        V++;
    }

    if(firstCodon == 'G' && secondCodon == 'U' && thirdCodon == 'G') {
        newGenome[0].codonUsage[36]++;
        V++;
    }

    if(firstCodon == 'G' && secondCodon == 'U' && thirdCodon == 'U') {
        newGenome[0].codonUsage[37]++;
        V++;
    }

    if(firstCodon == 'A' && secondCodon == 'A' && thirdCodon == 'A') {
        newGenome[0].codonUsage[38]++;
        K++;
    }

    if(firstCodon == 'A' && secondCodon == 'A' && thirdCodon == 'G') {
        newGenome[0].codonUsage[39]++;
        K++;
    }

    if(firstCodon == 'A' && secondCodon == 'A' && thirdCodon == 'C') {
        newGenome[0].codonUsage[40]++;
        N++;
    }

    if(firstCodon == 'A' && secondCodon == 'A' && thirdCodon == 'U') {
        newGenome[0].codonUsage[41]++;
        N++;
    }

    if(firstCodon == 'C' && secondCodon == 'A' && thirdCodon == 'A') {
        newGenome[0].codonUsage[42]++;
        Q++;
    }

    if(firstCodon == 'C' && secondCodon == 'A' && thirdCodon == 'G') {
        newGenome[0].codonUsage[43]++;
        Q++;
    }

    if(firstCodon == 'C' && secondCodon == 'A' && thirdCodon == 'C') {
        newGenome[0].codonUsage[44]++;
        H++;
    }

    if(firstCodon == 'C' && secondCodon == 'A' && thirdCodon == 'U') {
        newGenome[0].codonUsage[45]++;
        H++;
    }

    if(firstCodon == 'G' && secondCodon == 'A' && thirdCodon == 'A') {
        newGenome[0].codonUsage[46]++;
        E++;
    }

    if(firstCodon == 'G' && secondCodon == 'A' && thirdCodon == 'G') {
        newGenome[0].codonUsage[47]++;
        E++;
    }

    if(firstCodon == 'G' && secondCodon == 'A' && thirdCodon == 'C') {
        newGenome[0].codonUsage[48]++;
        D++;
    }

    if(firstCodon == 'G' && secondCodon == 'A' && thirdCodon == 'U') {
        newGenome[0].codonUsage[49]++;
        D++;
    }

    if(firstCodon == 'U' && secondCodon == 'A' && thirdCodon == 'C') {
        newGenome[0].codonUsage[50]++;
        Y++;
    }

    if(firstCodon == 'U' && secondCodon == 'A' && thirdCodon == 'U') {
        newGenome[0].codonUsage[51]++;
        Y++;
    }

    if(firstCodon == 'U' && secondCodon == 'G' && thirdCodon == 'C') {
        newGenome[0].codonUsage[52]++;
        C++;
    }

    if(firstCodon == 'U' && secondCodon == 'G' && thirdCodon == 'U') {
        newGenome[0].codonUsage[53]++;
        C++;
    }

    if(firstCodon == 'U' && secondCodon == 'U' && thirdCodon == 'C') {
        newGenome[0].codonUsage[54]++;
        F++;
    }

    if(firstCodon == 'U' && secondCodon == 'U' && thirdCodon == 'U') {
        newGenome[0].codonUsage[55]++;
        F++;
    }

    if(firstCodon == 'A' && secondCodon == 'U' && thirdCodon == 'A') {
        newGenome[0].codonUsage[56]++;
        I++;
    }

    if(firstCodon == 'A' && secondCodon == 'U' && thirdCodon == 'C') {
        newGenome[0].codonUsage[57]++;
        I++;
    }

    if(firstCodon == 'A' && secondCodon == 'U' && thirdCodon == 'U') {
        newGenome[0].codonUsage[58]++;
        I++;
    }

    if(firstCodon == 'A' && secondCodon == 'U' && thirdCodon == 'G') {
        newGenome[0].codonUsage[59]++;
        M++;
    }

    if(firstCodon == 'U' && secondCodon == 'G' && thirdCodon == 'G') {
        newGenome[0].codonUsage[60]++;
        W++;
    }

    if(firstCodon == 'U' && secondCodon == 'A' && thirdCodon == 'A') {
        newGenome[0].codonUsage[61]++;
        STOP++;
    }

    if(firstCodon == 'U' && secondCodon == 'A' && thirdCodon == 'G') {
        newGenome[0].codonUsage[62]++;
        STOP++;
    }

    if(firstCodon == 'U' && secondCodon == 'G' && thirdCodon == 'A') {
        newGenome[0].codonUsage[63]++;
        STOP++;
    }

// R,L,S,T,P,A,G,V,K,N,Q,H,E,D,Y,C,F,I,M,W,STOP
    inputProtein[0] = inputProtein[0] + R;
    inputProtein[1] = inputProtein[1] + L;
    inputProtein[2] = inputProtein[2] + S;
    inputProtein[3] = inputProtein[3] + T;
    inputProtein[4] = inputProtein[4] + P;
    inputProtein[5] = inputProtein[5] + A;
    inputProtein[6] = inputProtein[6] + G;
    inputProtein[7] = inputProtein[7] + V;
    inputProtein[8] = inputProtein[8] + K;
    inputProtein[9] = inputProtein[9] + N;
    inputProtein[10] = inputProtein[10] + Q;
    inputProtein[11] = inputProtein[11] + H;
    inputProtein[12] = inputProtein[12] + E;
    inputProtein[13] = inputProtein[13] + D;
    inputProtein[14] = inputProtein[14] + Y;
    inputProtein[15] = inputProtein[15] + C;
    inputProtein[16] = inputProtein[16] + F;
    inputProtein[17] = inputProtein[17] + I;
    inputProtein[18] = inputProtein[18] + M;
    inputProtein[19] = inputProtein[19] + W;
    inputProtein[20] = inputProtein[20] + STOP;


    return;
}


//*************************


void calculatePercentUsage(inputGenome newGenome[], int inputProtein[]) {


    for(int i = 0; i < 6; i++) {
        newGenome[0].proteinTotals[i] = (newGenome[0].codonUsage[i] / inputProtein[0]);
    }

    for(int i = 6; i < 12; i++) {
        newGenome[0].proteinTotals[i] = (newGenome[0].codonUsage[i] / inputProtein[1]);
    }

    for(int i = 12; i < 18; i++) {
        newGenome[0].proteinTotals[i] = (newGenome[0].codonUsage[i] / inputProtein[2]);
    }

    for(int i = 18; i < 22; i++) {
        newGenome[0].proteinTotals[i] = (newGenome[0].codonUsage[i] / inputProtein[3]);
    }

    for(int i = 22; i < 26; i++) {
        newGenome[0].proteinTotals[i] = (newGenome[0].codonUsage[i] / inputProtein[4]);
    }

    for(int i = 26; i < 30; i++) {
        newGenome[0].proteinTotals[i] = (newGenome[0].codonUsage[i] / inputProtein[5]);
    }

    for(int i = 30; i < 34; i++) {
        newGenome[0].proteinTotals[i] = (newGenome[0].codonUsage[i] / inputProtein[6]);
    }

    for(int i = 34; i < 38; i++) {
        newGenome[0].proteinTotals[i] = (newGenome[0].codonUsage[i] / inputProtein[7]);
    }

    for(int i = 38; i < 40; i++) {
        newGenome[0].proteinTotals[i] = (newGenome[0].codonUsage[i] / inputProtein[8]);
    }

    for(int i = 40; i < 42; i++) {
        newGenome[0].proteinTotals[i] = (newGenome[0].codonUsage[i] / inputProtein[9]);
    }

    for(int i = 42; i < 44; i++) {
        newGenome[0].proteinTotals[i] = (newGenome[0].codonUsage[i] / inputProtein[10]);
    }

    for(int i = 44; i < 46; i++) {
        newGenome[0].proteinTotals[i] = (newGenome[0].codonUsage[i] / inputProtein[11]);
    }

    for(int i = 46; i < 48; i++) {
        newGenome[0].proteinTotals[i] = (newGenome[0].codonUsage[i] / inputProtein[12]);
    }

    for(int i = 48; i < 50; i++) {
        newGenome[0].proteinTotals[i] = (newGenome[0].codonUsage[i] / inputProtein[13]);
    }

    for(int i = 50; i < 52; i++) {
        newGenome[0].proteinTotals[i] = (newGenome[0].codonUsage[i] / inputProtein[14]);
    }

    for(int i = 52; i < 54; i++) {
        newGenome[0].proteinTotals[i] = (newGenome[0].codonUsage[i] / inputProtein[15]);
    }

    for(int i = 54; i < 56; i++) {
        newGenome[0].proteinTotals[i] = (newGenome[0].codonUsage[i] / inputProtein[16]);
    }

    for(int i = 56; i < 59; i++) {
        newGenome[0].proteinTotals[i] = (newGenome[0].codonUsage[i] / inputProtein[17]);
    }

    newGenome[0].proteinTotals[59] = (newGenome[0].codonUsage[59] / inputProtein[18]);

    newGenome[0].proteinTotals[60] = (newGenome[0].codonUsage[60] / inputProtein[19]);

    for(int i = 61; i < 64; i++) {
        newGenome[0].proteinTotals[i] = (newGenome[0].codonUsage[i] / inputProtein[20]);
    }


    return;
}