#include <iostream>
#include <conio.h>
#include <ctime>
#include <cstdlib>
#include <windows.h>
#include <cctype>
#include <algorithm>

using namespace std;

int UInput;
const int columns = 11;
const int rows = 11;
char display[rows][columns];
char displaySpecial[rows][columns];
char displayPS[rows][columns];
char player = 'X';
char hit = '*';
char miss = 'M';
char tk;
int playerCol = 1;
int playerRow = 1;

char GetInput();


//Battleships.exe
//Introduction
void Introduction();
//(Give instructions and display options)

//Menu
void Menu();
void LoopBackToMenu();

//Option 1: Single Player
//Phase 1: Ship placement
//Phase 1 functions
void Phase1();
int RNGO();
void DisplayNormalGrid();

//Random Ship Placement (might overlap)
void DisplayScreenPS();
void RandomPlacement();

//global functions for the random placement function that will be used in phase 2
int PRNGR1(int AINumO1);
int PRNGC1(int AINumO1);

int PRNGR2(int AINumO2);
int PRNGC2(int AINumO2);

int PRNGR3(int AINumO3);
int PRNGC3(int AINumO3);

int PRNGR4(int AINumO4);
int PRNGC4(int AINumO4);

int PRNGR5(int AINumO5);
int PRNGC5(int AINumO5);

//Manual Ship Placement (might overlap)
void ManualPlacement();

//global variables from manual placement
int MNumR1;
int MNumC1;
int MNumO1;

int MNumR2;
int MNumC2;
int MNumO2;

int MNumR3;
int MNumC3;
int MNumO3;

int MNumR4;
int MNumC4;
int MNumO4;

int MNumR5;
int MNumC5;
int MNumO5;

//Phase 2
void Phase2(int PNumO1, int PNumR1, int PNumC1, int PNumO2, int PNumR2, int PNumC2, int PNumO3, int PNumR3, int PNumC3, int PNumO4, int PNumR4, int PNumC4, int PNumO5, int PNumR5, int PNumC5);
	//Phase 2 RNG number generator
	int AIRNGO1();
	int AIRNGR1(int AINumO1);
	int AIRNGC1(int AINumO1);

	int AIRNGO2();
	int AIRNGR2(int AINumO2);
	int AIRNGC2(int AINumO2);

	int AIRNGO3();
	int AIRNGR3(int AINumO3);
	int AIRNGC3(int AINumO3);

	int AIRNGO4();
	int AIRNGR4(int AINumO4);
	int AIRNGC4(int AINumO4);

	int AIRNGO5();
	int AIRNGR5(int AINumO5);
	int AIRNGC5(int AINumO5);

//The turn of the AI:
void AITurn(int AIShotR, int AIShotC, int PNumO1, int PNumR1, int PNumC1, int PNumO2, int PNumR2, int PNumC2, int PNumO3, int PNumR3, int PNumC3, int PNumO4, int PNumR4, int PNumC4, int PNumO5, int PNumR5, int PNumC5);
void PlayerTurn(int playerRow, int playerCol, int AINumO1, int AINumR1, int AINumC1, int AINumO2, int AINumR2, int AINumC2, int AINumO3, int AINumR3, int AINumC3, int AINumO4, int AINumR4, int AINumC4, int AINumO5, int AINumR5, int AINumC5);
//AI RNG for his turns
int AIRNGHR();
int AIRNGHC();


//Phase 2 grid display functions
void ClearScreen();
void DisplaySpecialGrid();
void DrawPlayer();
void DisplayScreen();

//Phase 3
//If player wins
void PlayerWins();
//If AI wins
void AIWins();
void Exit();

void Credits();

int main()
{
	system("cls");

	Introduction();

	Menu();
}

void Introduction()
{
	cout << "Welcome to Battleships.exe: " << endl;
	cout << "RULES: " << endl;
	cout << "The game has 3 phases." << endl;
	cout << "The 1st phase consists of placing your ships." << endl;
	cout << "You can place them manually one by one by placing in the coordinates" << endl;
	cout << "and orientation of the ship or you can get them randomly placed." << endl << endl;

	cout << "After you're done placing your ships, the 2nd phase begins in which" << endl;
	cout << "you get to shoot and get shot by the AI. In this phase, you move the" << endl;
	cout << "cursor and select where to shoot. You will be notified whether or not" << endl;
	cout << "you hit or missed an AI ship. After you shoot, the AI will randomly " << endl;
	cout << "shoot one of your grid tiles and notify you if he hit or missed." << endl << endl;

	cout << "Once the player or the AI has no more ships, the 3rd phase will begin;" << endl;
	cout << "where the winner will be determined." << endl << endl;

	cout << "To go into one of the following options, type in the number and then" << endl;
	cout << "press enter." << endl;
	cout << "This are the following options: " << endl;
	cout << "Option 1: Play singleplayer." << endl;
	cout << "Option 2: Exit the program." << endl;
	cout << "Option 3: See credits. " << endl;
	cout << "Enter your option: ";
}

void Menu()
{
	cin >> UInput;

	while (true)
	{
		if (UInput == 1)
		{
			Phase1();
		}
		else if (UInput == 2)
		{
			Exit();
		}
		else if (UInput == 3)
		{
			Credits();
		}
		else
		{
			cout << "Please enter a valid number: ";
			cin >> UInput;
		}
	}
}

void LoopBackToMenu()
{
	while (true)
	{
		if (UInput == 1)
		{
			main();
		}
		else
		{
			while (UInput != 1)
			{
				cout << "Please enter a valid number: ";
				cin >> UInput;
			}
			main();
		}
		return;
	}
}

//Phase 1
void Phase1() {
	srand(6);

	system("cls");
	cout << "This is Phase 1: SHIP PLACEMENT" << endl;
	cout << "Press 'R' or 'r' to place your ships randomly." << endl;
	cout << "Press 'M' or 'm' to place your ships manually." << endl;
	cout << "Would you like to place your ships manually or randomly?" << endl;
	cout << "Please enter one of those two letter: ";
	char tk = 0;
	while (tk != 'r' && tk != 'm')
	{
		tk = GetInput();
	}
	if (tk == 'r' || tk == 'R')
	{
		RandomPlacement();
	}
	if (tk == 'm' || tk == 'M')
	{
		ManualPlacement();
	}

	system("pause");
}

char GetInput() {
	while (_kbhit() == 0) {}
	char k = static_cast<char>(_getch());
	if (k > 0 && isalpha(k))
	{
		k = tolower(k);
	}
	return k;
}

int RNGO()
{
	int RNum2 = rand() % (2);
	return RNum2;
}

void DisplayNormalGrid()
{
	display[0][1] = '0'; display[0][2] = '1'; display[0][3] = '2'; display[0][4] = '3'; display[0][5] = '4';
	display[0][6] = '5'; display[0][7] = '6'; display[0][8] = '7'; display[0][9] = '8'; display[0][10] = '9';

	display[1][0] = '0'; display[2][0] = '1'; display[3][0] = '2'; display[4][0] = '3'; display[5][0] = '4';
	display[6][0] = '5'; display[7][0] = '6'; display[8][0] = '7'; display[9][0] = '8'; display[10][0] = '9';
}

void DisplayNormalGridPS()
{
	displayPS[0][1] = '0'; displayPS[0][2] = '1'; displayPS[0][3] = '2'; displayPS[0][4] = '3'; displayPS[0][5] = '4';
	displayPS[0][6] = '5'; displayPS[0][7] = '6'; displayPS[0][8] = '7'; displayPS[0][9] = '8'; displayPS[0][10] = '9';

	displayPS[1][0] = '0'; displayPS[2][0] = '1'; displayPS[3][0] = '2'; displayPS[4][0] = '3'; displayPS[5][0] = '4';
	displayPS[6][0] = '5'; displayPS[7][0] = '6'; displayPS[8][0] = '7'; displayPS[9][0] = '8'; displayPS[10][0] = '9';
}

void DisplayScreenPS() {
	for (int r = 0; r < rows; ++r) {
		for (int c = 0; c < columns; ++c) {
			std::cout << displayPS[r][c];
		}
		std::cout << std::endl;
	}
}

void RandomPlacement()
{
	system("cls");
	cout << "You chose to randomly place your ships." << endl;
	Sleep(500);

	//PB
	int PNumO1 = RNGO();
	int PNumR1 = PRNGR1(PNumO1);
	int PNumC1 = PRNGC1(PNumO1);

	int PNumO2 = RNGO();
	int PNumR2 = PRNGR2(PNumO2);
	int PNumC2 = PRNGC2(PNumO2);

	int PNumO3 = RNGO();
	int PNumR3 = PRNGR3(PNumO3);
	int PNumC3 = PRNGC3(PNumO3);

	int PNumO4 = RNGO();
	int PNumR4 = PRNGR4(PNumO4);
	int PNumC4 = PRNGC4(PNumO4);

	int PNumO5 = RNGO();
	int PNumR5 = PRNGR5(PNumO5);
	int PNumC5 = PRNGC5(PNumO5);

	if (PNumO1 == 0)
	{
		displayPS[PNumR1][PNumC1] = 'P';
		displayPS[PNumR1][PNumC1 + 1] = 'B';
	}
	else if (PNumO1 == 1)
	{
		displayPS[PNumR1][PNumC1] = 'P';
		displayPS[PNumR1 + 1][PNumC1] = 'B';
	}

	if (PNumO2 == 0)
	{
		displayPS[PNumR2][PNumC2] = 'S';
		displayPS[PNumR2][PNumC2 + 1] = 'U';
		displayPS[PNumR2][PNumC2 + 2] = 'B';
	}
	else if (PNumO2 == 1)
	{
		displayPS[PNumR2][PNumC2] = 'S';
		displayPS[PNumR2 + 1][PNumC2] = 'U';
		displayPS[PNumR2 + 2][PNumC2] = 'B';
	}

	if (PNumO3 == 0)
	{
		displayPS[PNumR3][PNumC3] = 'D';
		displayPS[PNumR3][PNumC3 + 1] = 'S';
		displayPS[PNumR3][PNumC3 + 2] = 'T';
	}
	else if (PNumO3 == 1)
	{
		displayPS[PNumR3][PNumC3] = 'D';
		displayPS[PNumR3 + 1][PNumC3] = 'S';
		displayPS[PNumR3 + 2][PNumC3] = 'T';
	}

	if (PNumO4 == 0)
	{
		displayPS[PNumR4][PNumC4] = 'B';
		displayPS[PNumR4][PNumC4 + 1] = 'T';
		displayPS[PNumR4][PNumC4 + 2] = 'S';
		displayPS[PNumR4][PNumC4 + 3] = 'P';
	}
	else if (PNumO4 == 1)
	{
		displayPS[PNumR4][PNumC4] = 'B';
		displayPS[PNumR4 + 1][PNumC4] = 'T';
		displayPS[PNumR4 + 2][PNumC4] = 'S';
		displayPS[PNumR4 + 3][PNumC4] = 'P';
	}

	if (PNumO5 == 0)
	{
		displayPS[PNumR5][PNumC5] = 'A';
		displayPS[PNumR5][PNumC5 + 1] = 'C';
		displayPS[PNumR5][PNumC5 + 2] = 'C';
		displayPS[PNumR5][PNumC5 + 3] = 'A';
		displayPS[PNumR5][PNumC5 + 4] = 'R';
	}
	else if (PNumO5 == 1)
	{
		displayPS[PNumR5][PNumC5] = 'A';
		displayPS[PNumR5 + 1][PNumC5] = 'C';
		displayPS[PNumR5 + 2][PNumC5] = 'C';
		displayPS[PNumR5 + 3][PNumC5] = 'A';
		displayPS[PNumR5 + 4][PNumC5] = 'R';
	}

	cout << "Your Ships: " << endl;
	DisplayNormalGridPS();
	DisplayScreenPS();

	cout << "Press 'y' or 'Y' when you are ready:" << endl;
	cout << "Are you ready? ";
	char tk1 = 0;
	while (tk1 != 'y') {
		tk1 = GetInput();
	}

	if (tk1 == 'y') {
		Phase2(PNumO1, PNumR1, PNumC1, PNumO2, PNumR2, PNumC2, PNumO3, PNumR3, PNumC3, PNumO4, PNumR4, PNumC4, PNumO5, PNumR5, PNumC5);

		system("cls");
		cout << "You win!" << endl;
		cout << "Congratulations!" << endl;
	}

	system("pause");
}

int PRNGR1(int PNumO1)
{
	if (PNumO1 == 1)
	{
		int PNumR = 1 + rand() % (9);
		return (PNumR);
	}
	if (PNumO1 == 0)
	{
		int PNumR = 1 + rand() % (10);
		return (PNumR);
	}
}
int PRNGC1(int PNumO1)
{
	if (PNumO1 == 1)
	{
		int PNumR = 1 + rand() % (10);
		return (PNumR);
	}
	if (PNumO1 == 0)
	{
		int PNumR = 1 + rand() % (9);
		return (PNumR);
	}
}

int PRNGR2(int PNumO2)
{
	if (PNumO2 == 1)
	{
		int PNumR = 1 + rand() % (8);
		return (PNumR);
	}
	if (PNumO2 == 0)
	{
		int PNumR = 1 + rand() % (10);
		return (PNumR);
	}
}
int PRNGC2(int PNumO2)
{
	if (PNumO2 == 1)
	{
		int PNumC = 1 + rand() % (10);
		return (PNumC);
	}
	if (PNumO2 == 0)
	{
		int PNumC = 1 + rand() % (8);
		return (PNumC);
	}
}

int PRNGR3(int PNumO3)
{
	if (PNumO3 == 1)
	{
		int PNumR = 1 + rand() % (8);
		return (PNumR);
	}
	if (PNumO3 == 0)
	{
		int PNumR = 1 + rand() % (10);
		return (PNumR);
	}
}
int PRNGC3(int PNumO3)
{
	if (PNumO3 == 1)
	{
		int PNumC = 1 + rand() % (10);
		return (PNumC);
	}
	if (PNumO3 == 0)
	{
		int PNumC = 1 + rand() % (8);
		return (PNumC);
	}
}

int PRNGR4(int PNumO4)
{
	if (PNumO4 == 1)
	{
		int PNumR = 1 + rand() % (7);
		return (PNumR);
	}
	if (PNumO4 == 0)
	{
		int PNumR = 1 + rand() % (10);
		return (PNumR);
	}
}
int PRNGC4(int PNumO4)
{
	if (PNumO4 == 1)
	{
		int PNumC = 1 + rand() % (10);
		return (PNumC);
	}
	if (PNumO4 == 0)
	{
		int PNumC = 1 + rand() % (7);
		return (PNumC);
	}
}

int PRNGR5(int PNumO5)
{
	if (PNumO5 == 1)
	{
		int PNumR = 1 + rand() % (6);
		return (PNumR);
	}
	if (PNumO5 == 0)
	{
		int PNumR = 1 + rand() % (10);
		return (PNumR);
	}
}
int PRNGC5(int PNumO5)
{
	if (PNumO5 == 1)
	{
		int PNumC = 1 + rand() % (10);
		return (PNumC);
	}
	if (PNumO5 == 0)
	{
		int PNumC = 1 + rand() % (6);
		return (PNumC);
	}
}

void ManualPlacement()
{
	system("cls");
	cout << "You chose to manually place your ships." << endl;
	cout << "The row and column range is from 0 to 9. " << endl;

	//PB
	int MNumR01;
	int MNumC01;

	cout << "The PATROL BOAT is 2 spaces long." << endl;
	cout << "Enter the orientation of the ship, 1 being vertical and 0 being horizontal: ";
	cin >> MNumO1;

	cout << "Enter the row you want the PATROL BOAT in: ";
	cin >> MNumR01;

	int MNumR1 = MNumR01 + 1;

	cout << "Enter the column you want the PATROL BOAT in: ";
	cin >> MNumC01;

	int MNumC1 = MNumC01 + 1;

	if (MNumO1 == 0)
	{
		displayPS[MNumR1][MNumC1] = 'P';
		if (MNumO1 == 1)
		{
			displayPS[MNumR1 + 1][MNumC1] = 'B';
		}
		else if (MNumO1 == 0)
		{
			displayPS[MNumR1][MNumC1 + 1] = 'B';
		}
	}
	else if (MNumO1 == 1)
	{
		displayPS[MNumR1][MNumC1] = 'P';
		if (MNumO1 == 1)
		{
			displayPS[MNumR1 + 1][MNumC1] = 'B';
		}
		else if (MNumO1 == 0)
		{
			displayPS[MNumR1][MNumC1 + 1] = 'B';
		}
	}

	//SUB
	system("cls");
	int MNumR02;
	int MNumC02;

	cout << "The SUBMARINE is 3 spaces long." << endl;
	cout << "Enter the orientation of the ship, 1 being vertical and 0 being horizontal: ";
	cin >> MNumO2;

	cout << "Enter the row you want the SUBMARINE in: ";
	cin >> MNumR02;

	int MNumR2 = MNumR02 + 1;

	cout << "Enter the column you want the SUBMARINE in: ";
	cin >> MNumC02;

	int MNumC2 = MNumC02 + 1;

	if (MNumO2 == 0)
	{
		displayPS[MNumR2][MNumC2] = 'S';
		if (MNumO2 == 1)
		{
			displayPS[MNumR2 + 1][MNumC2] = 'U';
			displayPS[MNumR2 + 2][MNumC2] = 'B';
		}
		else if (MNumO2 == 0)
		{
			displayPS[MNumR2][MNumC2 + 1] = 'U';
			displayPS[MNumR2][MNumC2 + 2] = 'B';
		}
	}
	else if (MNumO2 == 1)
	{
		displayPS[MNumR2][MNumC2] = 'S';
		if (MNumO2 == 1)
		{
			displayPS[MNumR2 + 1][MNumC2] = 'U';
			displayPS[MNumR2 + 2][MNumC2] = 'B';
		}
		else if (MNumO2 == 0)
		{
			displayPS[MNumR2][MNumC2 + 1] = 'U';
			displayPS[MNumR2][MNumC2 + 2] = 'B';
		}
	}

	//DST
	system("cls");
	int MNumR03;
	int MNumC03;

	cout << "The DESTROYER is 3 spaces long." << endl;
	cout << "Enter the orientation of the ship, 1 being vertical and 0 being horizontal: ";
	cin >> MNumO3;

	cout << "Enter the row you want the DESTROYER in: ";
	cin >> MNumR03;

	int MNumR3 = MNumR03 + 1;

	cout << "Enter the column you want the DESTROYER in: ";
	cin >> MNumC03;

	int MNumC3 = MNumC03 + 1;

	if (MNumO3 == 0)
	{
		displayPS[MNumR3][MNumC3] = 'D';
		if (MNumO3 == 1)
		{
			displayPS[MNumR3 + 1][MNumC3] = 'S';
			displayPS[MNumR3 + 2][MNumC3] = 'T';
		}
		else if (MNumO3 == 0)
		{
			displayPS[MNumR3][MNumC3 + 1] = 'S';
			displayPS[MNumR3][MNumC3 + 2] = 'T';
		}
	}
	else if (MNumO3 == 1)
	{
		displayPS[MNumR3][MNumC3] = 'D';
		if (MNumO3 == 1)
		{
			displayPS[MNumR3 + 1][MNumC3] = 'S';
			displayPS[MNumR3 + 2][MNumC3] = 'T';
		}
		else if (MNumO3 == 0)
		{
			displayPS[MNumR3][MNumC3 + 1] = 'S';
			displayPS[MNumR3][MNumC3 + 2] = 'T';
		}
	}

	//BTSP
	system("cls");
	int MNumR04;
	int MNumC04;

	cout << "The BATTLESHIP is 4 spaces long." << endl;
	cout << "Enter the orientation of the ship, 1 being vertical and 0 being horizontal: ";
	cin >> MNumO4;

	cout << "Enter the row you want the BATTLESHIP in: ";
	cin >> MNumR04;

	int MNumR4 = MNumR04 + 1;

	cout << "Enter the column you want the BATTLESHIP in: ";
	cin >> MNumC04;

	int MNumC4 = MNumC04 + 1;

	if (MNumO4 == 0)
	{
		displayPS[MNumR4][MNumC4] = 'B';
		if (MNumO4 == 1)
		{
			displayPS[MNumR4 + 1][MNumC4] = 'T';
			displayPS[MNumR4 + 2][MNumC4] = 'S';
			displayPS[MNumR4 + 3][MNumC4] = 'P';
		}
		else if (MNumO4 == 0)
		{
			displayPS[MNumR4][MNumC4 + 1] = 'T';
			displayPS[MNumR4][MNumC4 + 2] = 'S';
			displayPS[MNumR4][MNumC4 + 3] = 'P';
		}
	}
	else if (MNumO4 == 1)
	{
		displayPS[MNumR4][MNumC4] = 'B';
		if (MNumO4 == 1)
		{
			displayPS[MNumR4 + 1][MNumC4] = 'T';
			displayPS[MNumR4 + 2][MNumC4] = 'S';
			displayPS[MNumR4 + 3][MNumC4] = 'P';
		}
		else if (MNumO4 == 0)
		{
			displayPS[MNumR4][MNumC4 + 1] = 'T';
			displayPS[MNumR4][MNumC4 + 2] = 'S';
			displayPS[MNumR4][MNumC4 + 3] = 'P';
		}
	}

	//ACCAR
	system("cls");
	int MNumR05;
	int MNumC05;

	cout << "The AIRCRAFT CARRIER is 5 spaces long." << endl;
	cout << "Enter the orientation of the ship, 1 being vertical and 0 being horizontal: ";
	cin >> MNumO5;

	cout << "Enter the row you want the AIRCRAFT CARRIER in: ";
	cin >> MNumR05;

	int MNumR5 = MNumR05 + 1;

	cout << "Enter the column you want the AIRCRAFT CARRIER in: ";
	cin >> MNumC05;

	int MNumC5 = MNumC05 + 1;

	cout << endl;

	if (MNumO5 == 0)
	{
		displayPS[MNumR5][MNumC5] = 'A';
		if (MNumO5 == 1)
		{
			displayPS[MNumR5 + 1][MNumC5] = 'C';
			displayPS[MNumR5 + 2][MNumC5] = 'C';
			displayPS[MNumR5 + 3][MNumC5] = 'A';
			displayPS[MNumR5 + 4][MNumC5] = 'R';
		}
		else if (MNumO5 == 0)
		{
			displayPS[MNumR5][MNumC5 + 1] = 'C';
			displayPS[MNumR5][MNumC5 + 2] = 'C';
			displayPS[MNumR5][MNumC5 + 3] = 'A';
			displayPS[MNumR5][MNumC5 + 4] = 'R';
		}
	}
	else if (MNumO5 == 1)
	{
		displayPS[MNumR5][MNumC5] = 'A';
		if (MNumO5 == 1)
		{
			displayPS[MNumR5 + 1][MNumC5] = 'C';
			displayPS[MNumR5 + 2][MNumC5] = 'C';
			displayPS[MNumR5 + 3][MNumC5] = 'A';
			displayPS[MNumR5 + 4][MNumC5] = 'R';
		}
		else if (MNumO5 == 0)
		{
			displayPS[MNumR5][MNumC5 + 1] = 'C';
			displayPS[MNumR5][MNumC5 + 2] = 'C';
			displayPS[MNumR5][MNumC5 + 3] = 'A';
			displayPS[MNumR5][MNumC5 + 4] = 'R';
		}
	}

	cout << "Your Ships: " << endl;

	DisplayNormalGridPS();
	DisplayScreenPS();

	cout << "Press 'y' or 'Y' when you are ready:" << endl;
	cout << "Are you ready? ";
	char tk1 = 0;
	while (tk1 != 'y') {
		tk1 = GetInput();
	}

	if (tk1 == 'y') {
		Phase2(MNumO1, MNumR1, MNumC1, MNumO2, MNumR2, MNumC2, MNumO3, MNumR3, MNumC3, MNumO4, MNumR4, MNumC4, MNumO5, MNumR5, MNumC5);

		system("cls");
		cout << "You win!" << endl;
		cout << "Congratulations!" << endl;
	}

	system("pause");
}

//Phase 2
void Phase2(int PNumO1, int PNumR1, int PNumC1, int PNumO2, int PNumR2, int PNumC2, int PNumO3, int PNumR3, int PNumC3, int PNumO4, int PNumR4, int PNumC4, int PNumO5, int PNumR5, int PNumC5)
{
	memset(displaySpecial, ' ', sizeof(char) * columns * rows);
	bool bQuit = false;

	bool bSPB = false;
	bool bSPBTXT = true;

	bool bSSUB = false;
	bool bSSUBTXT = true;

	bool bSDST = false;
	bool bSDSTTXT = true;

	bool bSBTSP = false;
	bool bSBTSPTXT = true;

	bool bSACCAR = false;
	bool bSACCARTXT = true;

	int AINumO1 = AIRNGO1();
	int AINumR1 = AIRNGR1(AINumO1);
	int AINumC1 = AIRNGR1(AINumO1);

	int AINumO2 = AIRNGO2();
	int AINumR2 = AIRNGR2(AINumO2);
	int AINumC2 = AIRNGR2(AINumO2);

	int AINumO3 = AIRNGO3();
	int AINumR3 = AIRNGR3(AINumO3);
	int AINumC3 = AIRNGR3(AINumO3);

	int AINumO4 = AIRNGO4();
	int AINumR4 = AIRNGR4(AINumO4);
	int AINumC4 = AIRNGR4(AINumO4);

	int AINumO5 = AIRNGO5();
	int AINumR5 = AIRNGR5(AINumO5);
	int AINumC5 = AIRNGR5(AINumO5);

	bool bAISPB = false;
	bool bAISSUB = false;
	bool bAISDST = false;
	bool bAISBTSP = false;
	bool bAISACCAR = false;

	bool bSAIPBTXT = false;
	bool bSAISUBTXT = false;
	bool bSAIDSTTXT = false;
	bool bSAIBTSPTXT = false;
	bool bSAIACCARTXT = false;

	while (bQuit == false)
	{
		//Enemy Ship Grid:
		ClearScreen();
		DisplaySpecialGrid();
		DrawPlayer();
		DisplayScreen();

		//Text for hit and miss
		if (displaySpecial[playerRow][playerCol] == hit)
		{
			cout << "You hit an enemy ship!";
		}
		if (displaySpecial[playerRow][playerCol] == miss)
		{
			cout << "You missed!";
		}

		//PB Text
		if (AINumO1 == 1)
		{
			if ((displaySpecial[AINumR1][AINumC1] == hit) && (displaySpecial[AINumR1 + 1][AINumC1] == hit) && bSPBTXT)
			{
				cout << " and You sunk the enemy PATROL BOAT!";
				bSPBTXT = false;
			}
			else
			{
				cout << "";
			}
		}
		else if (AINumO1 == 0)
		{
			if ((displaySpecial[AINumR1][AINumC1] == hit) && (displaySpecial[AINumR1][AINumC1 + 1] == hit) && bSPBTXT)
			{
				cout << " and You sunk the enemy PATROL BOAT!";
				bSPBTXT = false;
			}
			else
			{
				cout << "";
			}
		}

		//SUB Text
		if (AINumO2 == 1)
		{
			if ((displaySpecial[AINumR2][AINumC2] == hit) && (displaySpecial[AINumR2 + 1][AINumC2] == hit) && (displaySpecial[AINumR2 + 2][AINumC2] == hit) && bSSUBTXT)
			{
				cout << " and You sunk the enemy SUBMARINE!";
				bSSUBTXT = false;
			}
			else
			{
				cout << "";
			}
		}
		else if (AINumO2 == 0)
		{
			if ((displaySpecial[AINumR2][AINumC2] == hit) && (displaySpecial[AINumR2][AINumC2 + 1] == hit) && (displaySpecial[AINumR2][AINumC2 + 2] == hit) && bSSUBTXT)
			{
				cout << " and You sunk the enemy SUBMARINE!";
				bSSUBTXT = false;
			}
			else
			{
				cout << "";
			}
		}

		//DST Text
		if (AINumO3 == 1)
		{
			if ((displaySpecial[AINumR3][AINumC3] == hit) && (displaySpecial[AINumR3 + 1][AINumC3] == hit) && (displaySpecial[AINumR3 + 2][AINumC3] == hit) && bSDSTTXT)
			{
				cout << " and You sunk the enemy DESTROYER!";
				bSDSTTXT = false;
			}
			else
			{
				cout << "";
			}
		}
		else if (AINumO3 == 0)
		{
			if ((displaySpecial[AINumR3][AINumC3] == hit) && (displaySpecial[AINumR3][AINumC3 + 1] == hit) && (displaySpecial[AINumR3][AINumC3 + 2] == hit) && bSDSTTXT)
			{
				cout << " and You sunk the enemy DESTROYER!";
				bSDSTTXT = false;
			}
			else
			{
				cout << "";
			}
		}

		//BTSP Text
		if (AINumO4 == 1)
		{
			if ((displaySpecial[AINumR4][AINumC4] == hit) && (displaySpecial[AINumR4 + 1][AINumC4] == hit) && (displaySpecial[AINumR4 + 2][AINumC4] == hit) && (displaySpecial[AINumR4 + 3][AINumC4] == hit) && bSBTSPTXT)
			{
				cout << " and You sunk the enemy BATTLESHIP!";
				bSBTSPTXT = false;
			}
			else
			{
				cout << "";
			}
		}
		else if (AINumO4 == 0)
		{
			if ((displaySpecial[AINumR4][AINumC4] == hit) && (displaySpecial[AINumR4][AINumC4 + 1] == hit) && (displaySpecial[AINumR4][AINumC4 + 2] == hit) && (displaySpecial[AINumR4][AINumC4 + 3] == hit) && bSBTSPTXT)
			{
				cout << " and You sunk the enemy BATTLESHIP!";
				bSBTSPTXT = false;
			}
			else
			{
				cout << "";
			}
		}

		//ACCAR Text
		if (AINumO5 == 1)
		{
			if ((displaySpecial[AINumR5][AINumC5] == hit) && (displaySpecial[AINumR5 + 1][AINumC5] == hit) && (displaySpecial[AINumR5 + 2][AINumC5] == hit) && (displaySpecial[AINumR5 + 3][AINumC5] == hit) && (displaySpecial[AINumR5 + 4][AINumC5] == hit) && bSACCARTXT)
			{
				cout << " and You sunk the enemy AIRCRAFT CARRIER!";
				bSACCARTXT = false;
			}
			else
			{
				cout << "";
			}
		}
		else if (AINumO5 == 0)
		{
			if ((displaySpecial[AINumR5][AINumC5] == hit) && (displaySpecial[AINumR5][AINumC5 + 1] == hit) && (displaySpecial[AINumR5][AINumC5 + 2] == hit) && (displaySpecial[AINumR5][AINumC5 + 3] == hit) && (displaySpecial[AINumR5][AINumC5 + 4] == hit) && bSACCARTXT)
			{
				cout << " and You sunk the enemy AIRCRAFT CARRIER!";
				bSACCARTXT = false;
			}
			else
			{
				cout << "";
			}
		}

		//Condition to see if PB sank
		if (AINumO1 == 1)
		{
			if ((displaySpecial[AINumR1][AINumC1] == hit) && (displaySpecial[AINumR1 + 1][AINumC1] == hit))
			{
				bSPB = true;
			}
		}
		else if (AINumO1 == 0)
		{
			if ((displaySpecial[AINumR1][AINumC1] == hit) && (displaySpecial[AINumR1][AINumC1 + 1] == hit))
			{
				bSPB = true;
			}
		}

		//Condition to see if SUB sank
		if (AINumO2 == 1)
		{
			if ((displaySpecial[AINumR2][AINumC2] == hit) && (displaySpecial[AINumR2 + 1][AINumC2] == hit) && (displaySpecial[AINumR2 + 2][AINumC2] == hit))
			{
				bSSUB = true;
			}
		}
		else if (AINumO2 == 0)
		{
			if ((displaySpecial[AINumR2][AINumC2] == hit) && (displaySpecial[AINumR2][AINumC2 + 1] == hit) && (displaySpecial[AINumR2][AINumC2 + 2] == hit))
			{
				bSSUB = true;
			}
		}

		//Condition to see if DST sank
		if (AINumO3 == 1)
		{
			if ((displaySpecial[AINumR3][AINumC3] == hit) && (displaySpecial[AINumR3 + 1][AINumC3] == hit) && (displaySpecial[AINumR3 + 2][AINumC3] == hit))
			{
				bSDST = true;
			}
		}
		else if (AINumO3 == 0)
		{
			if ((displaySpecial[AINumR3][AINumC3] == hit) && (displaySpecial[AINumR3][AINumC3 + 1] == hit) && (displaySpecial[AINumR3][AINumC3 + 2] == hit))
			{
				bSDST = true;
			}
		}

		//Condition to see if BTSP sank
		if (AINumO4 == 1)
		{
			if ((displaySpecial[AINumR4][AINumC4] == hit) && (displaySpecial[AINumR4 + 1][AINumC4] == hit) && (displaySpecial[AINumR4 + 2][AINumC4] == hit) && (displaySpecial[AINumR4 + 3][AINumC4] == hit))
			{
				bSBTSP = true;
			}
		}
		else if (AINumO4 == 0)
		{
			if ((displaySpecial[AINumR4][AINumC4] == hit) && (displaySpecial[AINumR4][AINumC4 + 1] == hit) && (displaySpecial[AINumR4][AINumC4 + 2] == hit) && (displaySpecial[AINumR4][AINumC4 + 3] == hit))
			{
				bSBTSP = true;
			}
		}

		//Condition to see if ACCAR sank
		if (AINumO5 == 1)
		{
			if ((displaySpecial[AINumR5][AINumC5] == hit) && (displaySpecial[AINumR5 + 1][AINumC5] == hit) && (displaySpecial[AINumR5 + 2][AINumC5] == hit) && (displaySpecial[AINumR5 + 3][AINumC5] == hit) && (displaySpecial[AINumR5 + 4][AINumC5] == hit))
			{
				bSACCAR = true;
			}
		}
		else if (AINumO5 == 0)
		{
			if ((displaySpecial[AINumR5][AINumC5] == hit) && (displaySpecial[AINumR5][AINumC5 + 1] == hit) && (displaySpecial[AINumR5][AINumC5 + 2] == hit) && (displaySpecial[AINumR5][AINumC5 + 3] == hit) && (displaySpecial[AINumR5][AINumC5 + 4] == hit))
			{
				bSACCAR = true;
			}
		}

		if ((bSPB) && (bSSUB) && (bSDST) && (bSBTSP) && (bSACCAR))
		{
			PlayerWins();
		}

		int AIShotR = AIRNGHR();
		int AIShotC = AIRNGHC();

		//Player's grid
		cout << endl;
		cout << "Your Ships:" << endl;

		DisplayNormalGridPS();
		DisplayScreenPS();

		char tk = GetInput();
		switch (tk)
		{
		case 'w':
			playerRow -= 1;
			if (playerRow == 0) {
				playerRow += 1;
			}
			break;
		case 's':
			playerRow += 1;
			if (playerRow > rows - 1) {
				playerRow -= 1;
			}
			break;
		case 'a':
			playerCol -= 1;
			if (playerCol == 0) {
				playerCol += 1;
			}
			break;
		case 'd':
			playerCol += 1;
			if (playerCol > columns - 1) {
				playerCol -= 1;
			}
			break;
		case 13:

			AITurn(AIShotR, AIShotC, PNumO1, PNumR1, PNumC1, PNumO2, PNumR2, PNumC2, PNumO3, PNumR3, PNumC3, PNumO4, PNumR4, PNumC4, PNumO5, PNumR5, PNumC5);



			if (PNumO1 == 1)
			{
				if ((displayPS[PNumR1][PNumC1] == hit) && (displayPS[PNumR1 + 1][PNumC1] == hit))
				{
					bAISPB = true;
				}
			}
			else if (PNumO1 == 0)
			{
				if ((displayPS[PNumR1][PNumC1] == hit) && (displayPS[PNumR1][PNumC1 + 1] == hit))
				{
					bAISPB = true;
				}
			}

			//Condition to see if SUbAI sank
			if (PNumO2 == 1)
			{
				if ((displayPS[PNumR2][PNumC2] == hit) && (displayPS[PNumR2 + 1][PNumC2] == hit) && (displayPS[PNumR2 + 2][PNumC2] == hit))
				{
					bAISSUB = true;
				}
			}
			else if (PNumO2 == 0)
			{
				if ((displayPS[PNumR2][PNumC2] == hit) && (displayPS[PNumR2][PNumC2 + 1] == hit) && (displayPS[PNumR2][PNumC2 + 2] == hit))
				{
					bAISSUB = true;
				}
			}

			//Condition to see if DST sank
			if (PNumO3 == 1)
			{
				if ((displayPS[PNumR3][PNumC3] == hit) && (displayPS[PNumR3 + 1][PNumC3] == hit) && (displayPS[PNumR3 + 2][PNumC3] == hit))
				{
					bAISDST = true;
				}
			}
			else if (PNumO3 == 0)
			{
				if ((displayPS[PNumR3][PNumC3] == hit) && (displayPS[PNumR3][PNumC3 + 1] == hit) && (displayPS[PNumR3][PNumC3 + 2] == hit))
				{
					bAISDST = true;
				}
			}

			//Condition to see if bAITSP sank
			if (PNumO4 == 1)
			{
				if ((displayPS[PNumR4][PNumC4] == hit) && (displayPS[PNumR4 + 1][PNumC4] == hit) && (displayPS[PNumR4 + 2][PNumC4] == hit) && (displayPS[PNumR4 + 3][PNumC4] == hit))
				{
					bAISBTSP = true;
				}
			}
			else if (PNumO4 == 0)
			{
				if ((displayPS[PNumR4][PNumC4] == hit) && (displayPS[PNumR4][PNumC4 + 1] == hit) && (displayPS[PNumR4][PNumC4 + 2] == hit) && (displayPS[PNumR4][PNumC4 + 3] == hit))
				{
					bAISBTSP = true;
				}
			}

			//Condition to see if ACCAR sank
			if (PNumO5 == 1)
			{
				if ((displayPS[PNumR5][PNumC5] == hit) && (displayPS[PNumR5 + 1][PNumC5] == hit) && (displayPS[PNumR5 + 2][PNumC5] == hit) && (displayPS[PNumR5 + 3][PNumC5] == hit) && (displayPS[PNumR5 + 4][PNumC5] == hit))
				{
					bAISACCAR = true;
				}
			}
			else if (PNumO5 == 0)
			{
				if ((displayPS[PNumR5][PNumC5] == hit) && (displayPS[PNumR5][PNumC5 + 1] == hit) && (displayPS[PNumR5][PNumC5 + 2] == hit) && (displayPS[PNumR5][PNumC5 + 3] == hit) && (displayPS[PNumR5][PNumC5 + 4] == hit))
				{
					bAISACCAR = true;
				}
			}

			if ((bAISPB) && (bAISSUB) && (bAISDST) && (bAISBTSP) && (bAISACCAR))
			{
				AIWins();
			}

			//Text for hit and miss for the AI
			if (displayPS[AIShotR][AIShotC] == hit)
			{
				cout << "The enemy has hit one of your ships!";
			}
			else
			{
				cout << "The enemy has missed one of your ships!";
			}
			Sleep(500);

			//All possibilites where the ships could be (might overlap)

			PlayerTurn(playerRow, playerCol, AINumO1, AINumR1, AINumC1, AINumO2, AINumR2, AINumC2, AINumO3, AINumR3, AINumC3, AINumO4, AINumR4, AINumC4, AINumO5, AINumR5, AINumC5);

			break;
		case 'q':
		case 27:
			bQuit = true;
			break;
		}
	}
}

//AIRNG1
int AIRNGO1()
{
	int AINumO1 = rand() % (2);
	return (AINumO1);
}
int AIRNGR1(int AINumO1)
{
	if (AINumO1 == 0)
	{
		int AINumR1 = 1 + rand() % (9);
		return (AINumR1);
	}
	if (AINumO1 == 1)
	{
		int AINumR1 = 1 + rand() % (10);
		return (AINumR1);
	}
}
int AIRNGC1(int AINumO1)
{
	if (AINumO1 == 0)
	{
		int AINumC1 = 1 + rand() % (10);
		return (AINumC1);
	}
	if (AINumO1 == 1)
	{
		int AINumC1 = 1 + rand() % (9);
		return (AINumC1);
	}
}

//AIRNG2
int AIRNGO2()
{
	int AINumO1 = rand() % (2);
	return (AINumO1);
}
int AIRNGR2(int AINumO2)
{
	if (AINumO2 == 0)
	{
		int AINumR1 = 1 + rand() % (8);
		return (AINumR1);
	}
	else if (AINumO2 == 1)
	{
		int AINumR1 = 1 + rand() % (10);
		return (AINumR1);
	}
}
int AIRNGC2(int AINumO2)
{
	if (AINumO2 == 0)
	{
		int AINumC1 = 1 + rand() % (10);
		return (AINumC1);
	}
	else if (AINumO2 == 1)
	{
		int AINumC1 = 1 + rand() % (8);
		return (AINumC1);
	}
}

//AIRNG3
int AIRNGO3()
{
	int AINumO3 = rand() % (2);
	return (AINumO3);
}
int AIRNGR3(int AINumO3)
{
	if (AINumO3 == 0)
	{
		int AINumR3 = 1 + rand() % (8);
		return (AINumR3);
	}
	else if (AINumO3 == 1)
	{
		int AINumR3 = 1 + rand() % (10);
		return (AINumR3);
	}
}
int AIRNGC3(int AINumO3)
{
	if (AINumO3 == 0)
	{
		int AINumC3 = 1 + rand() % (10);
		return (AINumC3);
	}
	else if (AINumO3 == 1)
	{
		int AINumC3 = 1 + rand() % (8);
		return (AINumC3);
	}
}


//AIRNG4
int AIRNGO4()
{
	int AINumO4 = rand() % (2);
	return (AINumO4);
}
int AIRNGR4(int AINumO4)
{
	if (AINumO4 == 0)
	{
		int AINumR4 = 1 + rand() % (7);
		return (AINumR4);
	}
	else if (AINumO4 == 1)
	{
		int AINumR4 = 1 + rand() % (10);
		return (AINumR4);
	}
}
int AIRNGC4(int AINumO4)
{
	if (AINumO4 == 0)
	{
		int AINumC4 = 1 + rand() % (10);
		return (AINumC4);
	}
	else if (AINumO4 == 1)
	{
		int AINumC4 = 1 + rand() % (7);
		return (AINumC4);
	}
}

//AIRNG5
int AIRNGO5()
{
	int AINumO5 = rand() % (2);
	return (AINumO5);
}

int AIRNGR5(int AINumO5)
{
	if (AINumO5 == 0)
	{
		int AINumR5 = 1 + rand() % (6);
		return (AINumR5);
	}
	else if (AINumO5 == 1)
	{
		int AINumR5 = 1 + rand() % (10);
		return (AINumR5);
	}
}

int AIRNGC5(int AINumO5)
{
	if (AINumO5 == 0)
	{
		int AINumC5 = 1 + rand() % (10);
		return (AINumC5);
	}
	else if (AINumO5 == 1)
	{
		int AINumC5 = 1 + rand() % (6);
		return (AINumC5);
	}
}

void ClearScreen()
{
	system("cls");
	memset(display, ' ', sizeof(char) * columns * rows);

	std::cout << "Enemy ships:" << std::endl;

	// x axis
	display[0][1] = '0'; display[0][2] = '1'; display[0][3] = '2'; display[0][4] = '3'; display[0][5] = '4';
	display[0][6] = '5'; display[0][7] = '6'; display[0][8] = '7'; display[0][9] = '8'; display[0][10] = '9';
	// y axis
	display[1][0] = '0'; display[2][0] = '1'; display[3][0] = '2'; display[4][0] = '3'; display[5][0] = '4';
	display[6][0] = '5'; display[7][0] = '6'; display[8][0] = '7'; display[9][0] = '8'; display[10][0] = '9';
}

void DisplaySpecialGrid()
{
	for (int r = 0; r < rows; ++r) {
		for (int c = 0; c < columns; ++c) {
			if (displaySpecial[r][c] != ' ')
			{
				display[r][c] = displaySpecial[r][c];
			}
		}
	}
}

void DisplayScreen() {
	for (int r = 0; r < rows; ++r) {
		for (int c = 0; c < columns; ++c) {
			std::cout << display[r][c];
		}
		std::cout << std::endl;
	}
}

void DrawPlayer() {
	display[playerRow][playerCol] = player;
}

void AITurn(int AIShotR, int AIShotC, int PNumO1, int PNumR1, int PNumC1, int PNumO2, int PNumR2, int PNumC2, int PNumO3, int PNumR3, int PNumC3, int PNumO4, int PNumR4, int PNumC4, int PNumO5, int PNumR5, int PNumC5)
{
	if (PNumO1 == 1)
	{
		if (PNumO2 == 1)
		{
			if (PNumO3 == 1)
			{
				if (PNumO4 == 1)
				{
					if (PNumO5 == 1)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 + 1 && AIShotC == PNumC1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 1 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 2 && AIShotC == PNumC2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 1 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 2 && AIShotC == PNumC3)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 1 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 2 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 3 && AIShotC == PNumC4)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 1 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 2 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 3 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 4 && AIShotC == PNumC5)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
					else if (PNumO5 == 0)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 + 1 && AIShotC == PNumC1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 1 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 2 && AIShotC == PNumC2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 1 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 2 && AIShotC == PNumC3)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 1 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 2 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 3 && AIShotC == PNumC4)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 1) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 2) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 3) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 4)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
				}
				else if (PNumO4 == 0)
				{
					if (PNumO5 == 1)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 + 1 && AIShotC == PNumC1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 1 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 2 && AIShotC == PNumC2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 1 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 2 && AIShotC == PNumC3)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 1) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 2) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 3)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 1 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 2 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 3 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 4 && AIShotC == PNumC5)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
					else if (PNumO5 == 0)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 + 1 && AIShotC == PNumC1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 1 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 2 && AIShotC == PNumC2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 1 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 2 && AIShotC == PNumC3)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 1) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 2) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 3)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 1) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 2) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 3) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 4)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
				}
			}
			else if (PNumO3 == 0)
			{
				if (PNumO4 == 1)
				{
					if (PNumO5 == 1)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 + 1 && AIShotC == PNumC1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 1 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 2 && AIShotC == PNumC2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 1) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 2)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 1 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 2 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 3 && AIShotC == PNumC4)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 1 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 2 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 3 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 4 && AIShotC == PNumC5)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
					else if (PNumO5 == 0)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 + 1 && AIShotC == PNumC1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 1 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 2 && AIShotC == PNumC2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 1) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 2)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 1 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 2 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 3 && AIShotC == PNumC4)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 1) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 2) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 3) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 4)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
				}
				else if (PNumO4 == 0)
				{
					if (PNumO5 == 1)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 + 1 && AIShotC == PNumC1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 1 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 2 && AIShotC == PNumC2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 1) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 2)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 1) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 2) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 3)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 1 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 2 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 3 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 4 && AIShotC == PNumC5)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
					else if (PNumO5 == 0)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 + 1 && AIShotC == PNumC1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 1 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 2 && AIShotC == PNumC2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 1) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 2)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 1) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 2) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 3)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 1) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 2) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 3) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 4)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
				}
			}
		}
		else if (PNumO2 == 0)
		{
			if (PNumO3 == 1)
			{
				if (PNumO4 == 0)
				{
					if (PNumO5 == 1)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 + 1 && AIShotC == PNumC1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 1) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 1 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 2 && AIShotC == PNumC3)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 1) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 2) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 3)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 1 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 2 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 3 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 4 && AIShotC == PNumC5)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
					else if (PNumO5 == 0)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 + 1 && AIShotC == PNumC1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 1) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 1 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 2 && AIShotC == PNumC3)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 1) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 2) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 3)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 1) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 2) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 3) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 4)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
				}
				else if (PNumO4 == 1)
				{
					if (PNumO5 == 1)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 + 1 && AIShotC == PNumC1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 1) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 1 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 2 && AIShotC == PNumC3)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 1 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 2 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 3 && AIShotC == PNumC4)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 1 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 2 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 3 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 4 && AIShotC == PNumC5)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
					else if (PNumO5 == 0)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 + 1 && AIShotC == PNumC1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 1) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 1 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 2 && AIShotC == PNumC3)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 1 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 2 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 3 && AIShotC == PNumC4)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 1) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 2) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 3) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 4)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
				}
			}
			else if (PNumO3 == 0)
			{
				if (PNumO4 == 1)
				{
					if (PNumO5 == 1)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 + 1 && AIShotC == PNumC1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 1) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 1) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 2)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 1 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 2 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 3 && AIShotC == PNumC4)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 1 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 2 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 3 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 4 && AIShotC == PNumC5)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
					else if (PNumO5 == 0)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 + 1 && AIShotC == PNumC1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 1) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 1) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 2)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 1 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 2 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 3 && AIShotC == PNumC4)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 1) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 2) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 3) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 4)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
				}
				else if (PNumO4 == 0)
				{
					if (PNumO5 == 1)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 + 1 && AIShotC == PNumC1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 1) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 1) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 2)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 1) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 2) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 3)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 1 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 2 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 3 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 4 && AIShotC == PNumC5)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
					else if (PNumO5 == 0)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 + 1 && AIShotC == PNumC1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 1) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 1) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 2)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 1) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 2) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 3)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 1) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 2) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 3) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 4)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
				}
			}
		}
	}

	if (PNumO1 == 0)
	{
		if (PNumO2 == 0)
		{
			if (PNumO3 == 0)
			{
				if (PNumO4 == 1)
				{
					if (PNumO5 == 1)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 && AIShotC == PNumC1 + 1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 1) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 1) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 2)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 1 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 2 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 3 && AIShotC == PNumC4)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 1 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 2 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 3 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 4 && AIShotC == PNumC5)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
					else if (PNumO5 == 0)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 && AIShotC == PNumC1 + 1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 1) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 1) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 2)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 1 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 2 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 3 && AIShotC == PNumC4)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 1) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 2) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 3) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 4)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
				}
				else if (PNumO4 == 0)
				{
					if (PNumO5 == 1)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 && AIShotC == PNumC1 + 1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 1) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 1) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 2)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 1) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 2) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 3)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 1 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 2 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 3 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 4 && AIShotC == PNumC5)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
					else if (PNumO5 == 0)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 && AIShotC == PNumC1 + 1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 1) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 1) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 2)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 1) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 2) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 3)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 1) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 2) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 3) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 4)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
				}
			}
			else if (PNumO3 == 1)
			{
				if (PNumO4 == 1)
				{
					if (PNumO5 == 1)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 && AIShotC == PNumC1 + 1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 1) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 1 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 2 && AIShotC == PNumC3)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 1 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 2 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 3 && AIShotC == PNumC4)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 1 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 2 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 3 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 4 && AIShotC == PNumC5)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
					else if (PNumO5 == 0)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 && AIShotC == PNumC1 + 1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 1) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 1 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 2 && AIShotC == PNumC3)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 1 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 2 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 3 && AIShotC == PNumC4)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 1) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 2) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 3) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 4)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
				}
				else if (PNumO4 == 0)
				{
					if (PNumO5 == 1)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 && AIShotC == PNumC1 + 1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 1) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 1 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 2 && AIShotC == PNumC3)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 1) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 2) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 3)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 1 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 2 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 3 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 4 && AIShotC == PNumC5)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
					else if (PNumO5 == 0)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 && AIShotC == PNumC1 + 1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 1) || (AIShotR == PNumR2 && AIShotC == PNumC2 + 2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 1 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 2 && AIShotC == PNumC3)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 1) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 2) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 3)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 1) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 2) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 3) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 4)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
				}
			}
		}
		else if (PNumO2 == 1)
		{
			if (PNumO3 == 0)
			{
				if (PNumO4 == 1)
				{
					if (PNumO5 == 1)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 && AIShotC == PNumC1 + 1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 1 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 2 && AIShotC == PNumC2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 1) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 2)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 1 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 2 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 3 && AIShotC == PNumC4)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 1 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 2 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 3 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 4 && AIShotC == PNumC5)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
					else if (PNumO5 == 0)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 && AIShotC == PNumC1 + 1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 1 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 2 && AIShotC == PNumC2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 1) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 2)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 1 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 2 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 3 && AIShotC == PNumC4)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 1) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 2) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 3) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 4)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
				}
				else if (PNumO4 == 0)
				{
					if (PNumO5 == 1)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 && AIShotC == PNumC1 + 1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 1 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 2 && AIShotC == PNumC2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 1) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 2)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 1) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 2) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 3)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 1 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 2 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 3 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 4 && AIShotC == PNumC5)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
					else if (PNumO5 == 0)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 && AIShotC == PNumC1 + 1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 1 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 2 && AIShotC == PNumC2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 1) || (AIShotR == PNumR3 && AIShotC == PNumC3 + 2)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 1) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 2) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 3)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 1) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 2) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 3) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 4)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
				}
			}
			else if (PNumO3 == 1)
			{
				if (PNumO4 == 1)
				{
					if (PNumO5 == 1)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 && AIShotC == PNumC1 + 1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 1 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 2 && AIShotC == PNumC2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 2 && AIShotC == PNumC3)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 1 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 2 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 3 && AIShotC == PNumC4)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 1 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 2 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 3 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 4 && AIShotC == PNumC5)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
					else if (PNumO5 == 0)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 && AIShotC == PNumC1 + 1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 1 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 2 && AIShotC == PNumC2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 2 && AIShotC == PNumC3)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 1 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 2 && AIShotC == PNumC4) || (AIShotR == PNumR4 + 3 && AIShotC == PNumC4)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 1) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 2) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 3) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 4)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
				}
				else if (PNumO4 == 0)
				{
					if (PNumO5 == 1)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 && AIShotC == PNumC1 + 1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 1 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 2 && AIShotC == PNumC2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 2 && AIShotC == PNumC3)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 1) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 2) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 3)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 1 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 2 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 3 && AIShotC == PNumC5) || (AIShotR == PNumR5 + 4 && AIShotC == PNumC5)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
					else if (PNumO5 == 0)
					{
						if (((AIShotR == PNumR1 && AIShotC == PNumC1) || (AIShotR == PNumR1 && AIShotC == PNumC1 + 1)) ||
							((AIShotR == PNumR2 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 1 && AIShotC == PNumC2) || (AIShotR == PNumR2 + 2 && AIShotC == PNumC2)) ||
							((AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 && AIShotC == PNumC3) || (AIShotR == PNumR3 + 2 && AIShotC == PNumC3)) ||
							((AIShotR == PNumR4 && AIShotC == PNumC4) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 1) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 2) || (AIShotR == PNumR4 && AIShotC == PNumC4 + 3)) ||
							((AIShotR == PNumR5 && AIShotC == PNumC5) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 1) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 2) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 3) || (AIShotR == PNumR5 && AIShotC == PNumC5 + 4)))
						{
							displayPS[AIShotR][AIShotC] = hit;
						}
						else
						{
							displayPS[AIShotR][AIShotC] = miss;
						}
					}
				}
			}
		}
	}
}

void PlayerTurn(int playerRow, int playerCol, int AINumO1, int AINumR1, int AINumC1, int AINumO2, int AINumR2, int AINumC2, int AINumO3, int AINumR3, int AINumC3, int AINumO4, int AINumR4, int AINumC4, int AINumO5, int AINumR5, int AINumC5)
{

	if (AINumO1 == 1)
	{
		if (AINumO2 == 1)
		{
			if (AINumO3 == 1)
			{
				if (AINumO4 == 1)
				{
					if (AINumO5 == 1)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 + 1 && playerCol == AINumC1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 + 1 && playerCol == AINumC2) || (playerRow == AINumR2 + 2 && playerCol == AINumC2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 + 1 && playerCol == AINumC3) || (playerRow == AINumR3 + 2 && playerCol == AINumC3)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 + 1 && playerCol == AINumC4) || (playerRow == AINumR4 + 2 && playerCol == AINumC4) || (playerRow == AINumR4 + 3 && playerCol == AINumC4)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 + 1 && playerCol == AINumC5) || (playerRow == AINumR5 + 2 && playerCol == AINumC5) || (playerRow == AINumR5 + 3 && playerCol == AINumC5) || (playerRow == AINumR5 + 4 && playerCol == AINumC5)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
					else if (AINumO5 == 0)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 + 1 && playerCol == AINumC1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 + 1 && playerCol == AINumC2) || (playerRow == AINumR2 + 2 && playerCol == AINumC2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 + 1 && playerCol == AINumC3) || (playerRow == AINumR3 + 2 && playerCol == AINumC3)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 + 1 && playerCol == AINumC4) || (playerRow == AINumR4 + 2 && playerCol == AINumC4) || (playerRow == AINumR4 + 3 && playerCol == AINumC4)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 && playerCol == AINumC5 + 1) || (playerRow == AINumR5 && playerCol == AINumC5 + 2) || (playerRow == AINumR5 && playerCol == AINumC5 + 3) || (playerRow == AINumR5 && playerCol == AINumC5 + 4)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
				}
				else if (AINumO4 == 0)
				{
					if (AINumO5 == 1)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 + 1 && playerCol == AINumC1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 + 1 && playerCol == AINumC2) || (playerRow == AINumR2 + 2 && playerCol == AINumC2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 + 1 && playerCol == AINumC3) || (playerRow == AINumR3 + 2 && playerCol == AINumC3)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 && playerCol == AINumC4 + 1) || (playerRow == AINumR4 && playerCol == AINumC4 + 2) || (playerRow == AINumR4 && playerCol == AINumC4 + 3)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 + 1 && playerCol == AINumC5) || (playerRow == AINumR5 + 2 && playerCol == AINumC5) || (playerRow == AINumR5 + 3 && playerCol == AINumC5) || (playerRow == AINumR5 + 4 && playerCol == AINumC5)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
					else if (AINumO5 == 0)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 + 1 && playerCol == AINumC1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 + 1 && playerCol == AINumC2) || (playerRow == AINumR2 + 2 && playerCol == AINumC2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 + 1 && playerCol == AINumC3) || (playerRow == AINumR3 + 2 && playerCol == AINumC3)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 && playerCol == AINumC4 + 1) || (playerRow == AINumR4 && playerCol == AINumC4 + 2) || (playerRow == AINumR4 && playerCol == AINumC4 + 3)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 && playerCol == AINumC5 + 1) || (playerRow == AINumR5 && playerCol == AINumC5 + 2) || (playerRow == AINumR5 && playerCol == AINumC5 + 3) || (playerRow == AINumR5 && playerCol == AINumC5 + 4)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
				}
			}
			else if (AINumO3 == 0)
			{
				if (AINumO4 == 1)
				{
					if (AINumO5 == 1)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 + 1 && playerCol == AINumC1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 + 1 && playerCol == AINumC2) || (playerRow == AINumR2 + 2 && playerCol == AINumC2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 && playerCol == AINumC3 + 1) || (playerRow == AINumR3 && playerCol == AINumC3 + 2)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 + 1 && playerCol == AINumC4) || (playerRow == AINumR4 + 2 && playerCol == AINumC4) || (playerRow == AINumR4 + 3 && playerCol == AINumC4)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 + 1 && playerCol == AINumC5) || (playerRow == AINumR5 + 2 && playerCol == AINumC5) || (playerRow == AINumR5 + 3 && playerCol == AINumC5) || (playerRow == AINumR5 + 4 && playerCol == AINumC5)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
					else if (AINumO5 == 0)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 + 1 && playerCol == AINumC1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 + 1 && playerCol == AINumC2) || (playerRow == AINumR2 + 2 && playerCol == AINumC2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 && playerCol == AINumC3 + 1) || (playerRow == AINumR3 && playerCol == AINumC3 + 2)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 + 1 && playerCol == AINumC4) || (playerRow == AINumR4 + 2 && playerCol == AINumC4) || (playerRow == AINumR4 + 3 && playerCol == AINumC4)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 && playerCol == AINumC5 + 1) || (playerRow == AINumR5 && playerCol == AINumC5 + 2) || (playerRow == AINumR5 && playerCol == AINumC5 + 3) || (playerRow == AINumR5 && playerCol == AINumC5 + 4)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
				}
				else if (AINumO4 == 0)
				{
					if (AINumO5 == 1)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 + 1 && playerCol == AINumC1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 + 1 && playerCol == AINumC2) || (playerRow == AINumR2 + 2 && playerCol == AINumC2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 && playerCol == AINumC3 + 1) || (playerRow == AINumR3 && playerCol == AINumC3 + 2)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 && playerCol == AINumC4 + 1) || (playerRow == AINumR4 && playerCol == AINumC4 + 2) || (playerRow == AINumR4 && playerCol == AINumC4 + 3)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 + 1 && playerCol == AINumC5) || (playerRow == AINumR5 + 2 && playerCol == AINumC5) || (playerRow == AINumR5 + 3 && playerCol == AINumC5) || (playerRow == AINumR5 + 4 && playerCol == AINumC5)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
					else if (AINumO5 == 0)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 + 1 && playerCol == AINumC1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 + 1 && playerCol == AINumC2) || (playerRow == AINumR2 + 2 && playerCol == AINumC2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 && playerCol == AINumC3 + 1) || (playerRow == AINumR3 && playerCol == AINumC3 + 2)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 && playerCol == AINumC4 + 1) || (playerRow == AINumR4 && playerCol == AINumC4 + 2) || (playerRow == AINumR4 && playerCol == AINumC4 + 3)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 && playerCol == AINumC5 + 1) || (playerRow == AINumR5 && playerCol == AINumC5 + 2) || (playerRow == AINumR5 && playerCol == AINumC5 + 3) || (playerRow == AINumR5 && playerCol == AINumC5 + 4)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
				}
			}
		}
		else if (AINumO2 == 0)
		{
			if (AINumO3 == 1)
			{
				if (AINumO4 == 0)
				{
					if (AINumO5 == 1)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 + 1 && playerCol == AINumC1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 && playerCol == AINumC2 + 1) || (playerRow == AINumR2 && playerCol == AINumC2 + 2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 + 1 && playerCol == AINumC3) || (playerRow == AINumR3 + 2 && playerCol == AINumC3)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 && playerCol == AINumC4 + 1) || (playerRow == AINumR4 && playerCol == AINumC4 + 2) || (playerRow == AINumR4 && playerCol == AINumC4 + 3)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 + 1 && playerCol == AINumC5) || (playerRow == AINumR5 + 2 && playerCol == AINumC5) || (playerRow == AINumR5 + 3 && playerCol == AINumC5) || (playerRow == AINumR5 + 4 && playerCol == AINumC5)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
					else if (AINumO5 == 0)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 + 1 && playerCol == AINumC1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 && playerCol == AINumC2 + 1) || (playerRow == AINumR2 && playerCol == AINumC2 + 2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 + 1 && playerCol == AINumC3) || (playerRow == AINumR3 + 2 && playerCol == AINumC3)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 && playerCol == AINumC4 + 1) || (playerRow == AINumR4 && playerCol == AINumC4 + 2) || (playerRow == AINumR4 && playerCol == AINumC4 + 3)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 && playerCol == AINumC5 + 1) || (playerRow == AINumR5 && playerCol == AINumC5 + 2) || (playerRow == AINumR5 && playerCol == AINumC5 + 3) || (playerRow == AINumR5 && playerCol == AINumC5 + 4)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
				}
				else if (AINumO4 == 1)
				{
					if (AINumO5 == 1)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 + 1 && playerCol == AINumC1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 && playerCol == AINumC2 + 1) || (playerRow == AINumR2 && playerCol == AINumC2 + 2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 + 1 && playerCol == AINumC3) || (playerRow == AINumR3 + 2 && playerCol == AINumC3)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 + 1 && playerCol == AINumC4) || (playerRow == AINumR4 + 2 && playerCol == AINumC4) || (playerRow == AINumR4 + 3 && playerCol == AINumC4)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 + 1 && playerCol == AINumC5) || (playerRow == AINumR5 + 2 && playerCol == AINumC5) || (playerRow == AINumR5 + 3 && playerCol == AINumC5) || (playerRow == AINumR5 + 4 && playerCol == AINumC5)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
					else if (AINumO5 == 0)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 + 1 && playerCol == AINumC1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 && playerCol == AINumC2 + 1) || (playerRow == AINumR2 && playerCol == AINumC2 + 2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 + 1 && playerCol == AINumC3) || (playerRow == AINumR3 + 2 && playerCol == AINumC3)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 + 1 && playerCol == AINumC4) || (playerRow == AINumR4 + 2 && playerCol == AINumC4) || (playerRow == AINumR4 + 3 && playerCol == AINumC4)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 && playerCol == AINumC5 + 1) || (playerRow == AINumR5 && playerCol == AINumC5 + 2) || (playerRow == AINumR5 && playerCol == AINumC5 + 3) || (playerRow == AINumR5 && playerCol == AINumC5 + 4)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
				}
			}
			else if (AINumO3 == 0)
			{
				if (AINumO4 == 1)
				{
					if (AINumO5 == 1)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 + 1 && playerCol == AINumC1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 && playerCol == AINumC2 + 1) || (playerRow == AINumR2 && playerCol == AINumC2 + 2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 && playerCol == AINumC3 + 1) || (playerRow == AINumR3 && playerCol == AINumC3 + 2)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 + 1 && playerCol == AINumC4) || (playerRow == AINumR4 + 2 && playerCol == AINumC4) || (playerRow == AINumR4 + 3 && playerCol == AINumC4)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 + 1 && playerCol == AINumC5) || (playerRow == AINumR5 + 2 && playerCol == AINumC5) || (playerRow == AINumR5 + 3 && playerCol == AINumC5) || (playerRow == AINumR5 + 4 && playerCol == AINumC5)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
					else if (AINumO5 == 0)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 + 1 && playerCol == AINumC1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 && playerCol == AINumC2 + 1) || (playerRow == AINumR2 && playerCol == AINumC2 + 2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 && playerCol == AINumC3 + 1) || (playerRow == AINumR3 && playerCol == AINumC3 + 2)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 + 1 && playerCol == AINumC4) || (playerRow == AINumR4 + 2 && playerCol == AINumC4) || (playerRow == AINumR4 + 3 && playerCol == AINumC4)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 && playerCol == AINumC5 + 1) || (playerRow == AINumR5 && playerCol == AINumC5 + 2) || (playerRow == AINumR5 && playerCol == AINumC5 + 3) || (playerRow == AINumR5 && playerCol == AINumC5 + 4)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
				}
				else if (AINumO4 == 0)
				{
					if (AINumO5 == 1)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 + 1 && playerCol == AINumC1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 && playerCol == AINumC2 + 1) || (playerRow == AINumR2 && playerCol == AINumC2 + 2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 && playerCol == AINumC3 + 1) || (playerRow == AINumR3 && playerCol == AINumC3 + 2)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 && playerCol == AINumC4 + 1) || (playerRow == AINumR4 && playerCol == AINumC4 + 2) || (playerRow == AINumR4 && playerCol == AINumC4 + 3)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 + 1 && playerCol == AINumC5) || (playerRow == AINumR5 + 2 && playerCol == AINumC5) || (playerRow == AINumR5 + 3 && playerCol == AINumC5) || (playerRow == AINumR5 + 4 && playerCol == AINumC5)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
					else if (AINumO5 == 0)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 + 1 && playerCol == AINumC1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 && playerCol == AINumC2 + 1) || (playerRow == AINumR2 && playerCol == AINumC2 + 2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 && playerCol == AINumC3 + 1) || (playerRow == AINumR3 && playerCol == AINumC3 + 2)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 && playerCol == AINumC4 + 1) || (playerRow == AINumR4 && playerCol == AINumC4 + 2) || (playerRow == AINumR4 && playerCol == AINumC4 + 3)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 && playerCol == AINumC5 + 1) || (playerRow == AINumR5 && playerCol == AINumC5 + 2) || (playerRow == AINumR5 && playerCol == AINumC5 + 3) || (playerRow == AINumR5 && playerCol == AINumC5 + 4)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
				}
			}
		}
	}

	if (AINumO1 == 0)
	{
		if (AINumO2 == 0)
		{
			if (AINumO3 == 0)
			{
				if (AINumO4 == 1)
				{
					if (AINumO5 == 1)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 && playerCol == AINumC1 + 1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 && playerCol == AINumC2 + 1) || (playerRow == AINumR2 && playerCol == AINumC2 + 2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 && playerCol == AINumC3 + 1) || (playerRow == AINumR3 && playerCol == AINumC3 + 2)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 + 1 && playerCol == AINumC4) || (playerRow == AINumR4 + 2 && playerCol == AINumC4) || (playerRow == AINumR4 + 3 && playerCol == AINumC4)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 + 1 && playerCol == AINumC5) || (playerRow == AINumR5 + 2 && playerCol == AINumC5) || (playerRow == AINumR5 + 3 && playerCol == AINumC5) || (playerRow == AINumR5 + 4 && playerCol == AINumC5)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
					else if (AINumO5 == 0)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 && playerCol == AINumC1 + 1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 && playerCol == AINumC2 + 1) || (playerRow == AINumR2 && playerCol == AINumC2 + 2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 && playerCol == AINumC3 + 1) || (playerRow == AINumR3 && playerCol == AINumC3 + 2)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 + 1 && playerCol == AINumC4) || (playerRow == AINumR4 + 2 && playerCol == AINumC4) || (playerRow == AINumR4 + 3 && playerCol == AINumC4)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 && playerCol == AINumC5 + 1) || (playerRow == AINumR5 && playerCol == AINumC5 + 2) || (playerRow == AINumR5 && playerCol == AINumC5 + 3) || (playerRow == AINumR5 && playerCol == AINumC5 + 4)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
				}
				else if (AINumO4 == 0)
				{
					if (AINumO5 == 1)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 && playerCol == AINumC1 + 1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 && playerCol == AINumC2 + 1) || (playerRow == AINumR2 && playerCol == AINumC2 + 2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 && playerCol == AINumC3 + 1) || (playerRow == AINumR3 && playerCol == AINumC3 + 2)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 && playerCol == AINumC4 + 1) || (playerRow == AINumR4 && playerCol == AINumC4 + 2) || (playerRow == AINumR4 && playerCol == AINumC4 + 3)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 + 1 && playerCol == AINumC5) || (playerRow == AINumR5 + 2 && playerCol == AINumC5) || (playerRow == AINumR5 + 3 && playerCol == AINumC5) || (playerRow == AINumR5 + 4 && playerCol == AINumC5)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
					else if (AINumO5 == 0)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 && playerCol == AINumC1 + 1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 && playerCol == AINumC2 + 1) || (playerRow == AINumR2 && playerCol == AINumC2 + 2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 && playerCol == AINumC3 + 1) || (playerRow == AINumR3 && playerCol == AINumC3 + 2)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 && playerCol == AINumC4 + 1) || (playerRow == AINumR4 && playerCol == AINumC4 + 2) || (playerRow == AINumR4 && playerCol == AINumC4 + 3)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 && playerCol == AINumC5 + 1) || (playerRow == AINumR5 && playerCol == AINumC5 + 2) || (playerRow == AINumR5 && playerCol == AINumC5 + 3) || (playerRow == AINumR5 && playerCol == AINumC5 + 4)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
				}
			}
			else if (AINumO3 == 1)
			{
				if (AINumO4 == 1)
				{
					if (AINumO5 == 1)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 && playerCol == AINumC1 + 1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 && playerCol == AINumC2 + 1) || (playerRow == AINumR2 && playerCol == AINumC2 + 2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 + 1 && playerCol == AINumC3) || (playerRow == AINumR3 + 2 && playerCol == AINumC3)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 + 1 && playerCol == AINumC4) || (playerRow == AINumR4 + 2 && playerCol == AINumC4) || (playerRow == AINumR4 + 3 && playerCol == AINumC4)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 + 1 && playerCol == AINumC5) || (playerRow == AINumR5 + 2 && playerCol == AINumC5) || (playerRow == AINumR5 + 3 && playerCol == AINumC5) || (playerRow == AINumR5 + 4 && playerCol == AINumC5)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
					else if (AINumO5 == 0)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 && playerCol == AINumC1 + 1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 && playerCol == AINumC2 + 1) || (playerRow == AINumR2 && playerCol == AINumC2 + 2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 + 1 && playerCol == AINumC3) || (playerRow == AINumR3 + 2 && playerCol == AINumC3)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 + 1 && playerCol == AINumC4) || (playerRow == AINumR4 + 2 && playerCol == AINumC4) || (playerRow == AINumR4 + 3 && playerCol == AINumC4)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 && playerCol == AINumC5 + 1) || (playerRow == AINumR5 && playerCol == AINumC5 + 2) || (playerRow == AINumR5 && playerCol == AINumC5 + 3) || (playerRow == AINumR5 && playerCol == AINumC5 + 4)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
				}
				else if (AINumO4 == 0)
				{
					if (AINumO5 == 1)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 && playerCol == AINumC1 + 1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 && playerCol == AINumC2 + 1) || (playerRow == AINumR2 && playerCol == AINumC2 + 2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 + 1 && playerCol == AINumC3) || (playerRow == AINumR3 + 2 && playerCol == AINumC3)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 && playerCol == AINumC4 + 1) || (playerRow == AINumR4 && playerCol == AINumC4 + 2) || (playerRow == AINumR4 && playerCol == AINumC4 + 3)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 + 1 && playerCol == AINumC5) || (playerRow == AINumR5 + 2 && playerCol == AINumC5) || (playerRow == AINumR5 + 3 && playerCol == AINumC5) || (playerRow == AINumR5 + 4 && playerCol == AINumC5)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
					else if (AINumO5 == 0)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 && playerCol == AINumC1 + 1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 && playerCol == AINumC2 + 1) || (playerRow == AINumR2 && playerCol == AINumC2 + 2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 + 1 && playerCol == AINumC3) || (playerRow == AINumR3 + 2 && playerCol == AINumC3)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 && playerCol == AINumC4 + 1) || (playerRow == AINumR4 && playerCol == AINumC4 + 2) || (playerRow == AINumR4 && playerCol == AINumC4 + 3)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 && playerCol == AINumC5 + 1) || (playerRow == AINumR5 && playerCol == AINumC5 + 2) || (playerRow == AINumR5 && playerCol == AINumC5 + 3) || (playerRow == AINumR5 && playerCol == AINumC5 + 4)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
				}
			}
		}
		else if (AINumO2 == 1)
		{
			if (AINumO3 == 0)
			{
				if (AINumO4 == 1)
				{
					if (AINumO5 == 1)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 && playerCol == AINumC1 + 1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 + 1 && playerCol == AINumC2) || (playerRow == AINumR2 + 2 && playerCol == AINumC2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 && playerCol == AINumC3 + 1) || (playerRow == AINumR3 && playerCol == AINumC3 + 2)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 + 1 && playerCol == AINumC4) || (playerRow == AINumR4 + 2 && playerCol == AINumC4) || (playerRow == AINumR4 + 3 && playerCol == AINumC4)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 + 1 && playerCol == AINumC5) || (playerRow == AINumR5 + 2 && playerCol == AINumC5) || (playerRow == AINumR5 + 3 && playerCol == AINumC5) || (playerRow == AINumR5 + 4 && playerCol == AINumC5)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
					else if (AINumO5 == 0)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 && playerCol == AINumC1 + 1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 + 1 && playerCol == AINumC2) || (playerRow == AINumR2 + 2 && playerCol == AINumC2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 && playerCol == AINumC3 + 1) || (playerRow == AINumR3 && playerCol == AINumC3 + 2)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 + 1 && playerCol == AINumC4) || (playerRow == AINumR4 + 2 && playerCol == AINumC4) || (playerRow == AINumR4 + 3 && playerCol == AINumC4)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 && playerCol == AINumC5 + 1) || (playerRow == AINumR5 && playerCol == AINumC5 + 2) || (playerRow == AINumR5 && playerCol == AINumC5 + 3) || (playerRow == AINumR5 && playerCol == AINumC5 + 4)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
				}
				else if (AINumO4 == 0)
				{
					if (AINumO5 == 1)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 && playerCol == AINumC1 + 1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 + 1 && playerCol == AINumC2) || (playerRow == AINumR2 + 2 && playerCol == AINumC2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 && playerCol == AINumC3 + 1) || (playerRow == AINumR3 && playerCol == AINumC3 + 2)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 && playerCol == AINumC4 + 1) || (playerRow == AINumR4 && playerCol == AINumC4 + 2) || (playerRow == AINumR4 && playerCol == AINumC4 + 3)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 + 1 && playerCol == AINumC5) || (playerRow == AINumR5 + 2 && playerCol == AINumC5) || (playerRow == AINumR5 + 3 && playerCol == AINumC5) || (playerRow == AINumR5 + 4 && playerCol == AINumC5)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
					else if (AINumO5 == 0)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 && playerCol == AINumC1 + 1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 + 1 && playerCol == AINumC2) || (playerRow == AINumR2 + 2 && playerCol == AINumC2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 && playerCol == AINumC3 + 1) || (playerRow == AINumR3 && playerCol == AINumC3 + 2)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 && playerCol == AINumC4 + 1) || (playerRow == AINumR4 && playerCol == AINumC4 + 2) || (playerRow == AINumR4 && playerCol == AINumC4 + 3)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 && playerCol == AINumC5 + 1) || (playerRow == AINumR5 && playerCol == AINumC5 + 2) || (playerRow == AINumR5 && playerCol == AINumC5 + 3) || (playerRow == AINumR5 && playerCol == AINumC5 + 4)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
				}
			}
			else if (AINumO3 == 1)
			{
				if (AINumO4 == 1)
				{
					if (AINumO5 == 1)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 && playerCol == AINumC1 + 1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 + 1 && playerCol == AINumC2) || (playerRow == AINumR2 + 2 && playerCol == AINumC2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 + 2 && playerCol == AINumC3)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 + 1 && playerCol == AINumC4) || (playerRow == AINumR4 + 2 && playerCol == AINumC4) || (playerRow == AINumR4 + 3 && playerCol == AINumC4)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 + 1 && playerCol == AINumC5) || (playerRow == AINumR5 + 2 && playerCol == AINumC5) || (playerRow == AINumR5 + 3 && playerCol == AINumC5) || (playerRow == AINumR5 + 4 && playerCol == AINumC5)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
					else if (AINumO5 == 0)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 && playerCol == AINumC1 + 1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 + 1 && playerCol == AINumC2) || (playerRow == AINumR2 + 2 && playerCol == AINumC2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 + 2 && playerCol == AINumC3)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 + 1 && playerCol == AINumC4) || (playerRow == AINumR4 + 2 && playerCol == AINumC4) || (playerRow == AINumR4 + 3 && playerCol == AINumC4)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 && playerCol == AINumC5 + 1) || (playerRow == AINumR5 && playerCol == AINumC5 + 2) || (playerRow == AINumR5 && playerCol == AINumC5 + 3) || (playerRow == AINumR5 && playerCol == AINumC5 + 4)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
				}
				else if (AINumO4 == 0)
				{
					if (AINumO5 == 1)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 && playerCol == AINumC1 + 1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 + 1 && playerCol == AINumC2) || (playerRow == AINumR2 + 2 && playerCol == AINumC2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 + 2 && playerCol == AINumC3)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 && playerCol == AINumC4 + 1) || (playerRow == AINumR4 && playerCol == AINumC4 + 2) || (playerRow == AINumR4 && playerCol == AINumC4 + 3)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 + 1 && playerCol == AINumC5) || (playerRow == AINumR5 + 2 && playerCol == AINumC5) || (playerRow == AINumR5 + 3 && playerCol == AINumC5) || (playerRow == AINumR5 + 4 && playerCol == AINumC5)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
					else if (AINumO5 == 0)
					{
						if (((playerRow == AINumR1 && playerCol == AINumC1) || (playerRow == AINumR1 && playerCol == AINumC1 + 1)) ||
							((playerRow == AINumR2 && playerCol == AINumC2) || (playerRow == AINumR2 + 1 && playerCol == AINumC2) || (playerRow == AINumR2 + 2 && playerCol == AINumC2)) ||
							((playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 && playerCol == AINumC3) || (playerRow == AINumR3 + 2 && playerCol == AINumC3)) ||
							((playerRow == AINumR4 && playerCol == AINumC4) || (playerRow == AINumR4 && playerCol == AINumC4 + 1) || (playerRow == AINumR4 && playerCol == AINumC4 + 2) || (playerRow == AINumR4 && playerCol == AINumC4 + 3)) ||
							((playerRow == AINumR5 && playerCol == AINumC5) || (playerRow == AINumR5 && playerCol == AINumC5 + 1) || (playerRow == AINumR5 && playerCol == AINumC5 + 2) || (playerRow == AINumR5 && playerCol == AINumC5 + 3) || (playerRow == AINumR5 && playerCol == AINumC5 + 4)))
						{
							displaySpecial[playerRow][playerCol] = hit;
						}
						else
						{
							displaySpecial[playerRow][playerCol] = miss;
						}
					}
				}
			}
		}
	}
}

int AIRNGHR()
{
	int AIHR = 1 + rand() % 10;
	return (AIHR);
}

int AIRNGHC()
{
	int AIHC = 1 + rand() % 10;
	return (AIHC);
}

void PlayerWins()
{
	system("cls");
	cout << "Congratulations!" << endl;
	cout << "You beat the AI by sinking it's entire fleet!" << endl;
	cout << "If you would like to leave enter 1: ";
	cin >> UInput;
	if (UInput == 1)
	{
		Exit();
	}
}

void AIWins()
{
	system("cls");
	cout << "YOU LOSE!" << endl;
	cout << "The AI beat you by sinking your entire fleet!" << endl;
	cout << "If you would like to leave enter 1: ";
	cin >> UInput;
	if (UInput == 1)
	{
		Exit();
	}
}

void Exit() { exit(0); }

void Credits()
{
	system("cls");
	cout << "Made by Daniel Opoka." << endl;
	cout << "If you would like to go back to the main menu enter 1: ";
	cin >> UInput;
	if (UInput == 1)
	{
		main();
	}
	else
	{
		LoopBackToMenu();
	}
}