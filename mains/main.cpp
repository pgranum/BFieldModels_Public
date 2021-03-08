#include <iostream>

#include "BFieldModels.h"
#include "Loop.h"
#include "mcDTubeSupportClass.h"
#include "mirrorModels.h"
#include "Shell.h"
#include "Tube.h"

int main(){
	double carP[3] = {0.};
	double BCar[3] = {0.};
	
	Loop loop = loop();
	const int n = 7;
	
	mcDonald(loop,n,carP,BCar);
	
}
