#include <ctype.h>
#include <stdio.h>
#include "utils.h"

unsigned long atoul(char *a)
{
	unsigned long ul = 0;

	while (isdigit(*a))
	{
		ul = ul*10 + *a-48;
		a++;
	}

	return ul;
}
