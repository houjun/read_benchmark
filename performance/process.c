#include <stdio.h>

void main(int argc, char *argv[])
{
    
    FILE *fp = fopen(argv[1], "r");
    if (fp == NULL) {
        printf("Error in opening a file: %s", argv[1]);
        return;
    }

    char tmp[1024];
    char tmp2[1024];
    double p;
    int i=1;
    while(fgets (tmp ,1024 ,fp) != NULL) {
        if(sscanf(tmp, "Total time: %lf %s", &p, tmp2) == 2) {
            printf("%lf", p);
            if(i%2==1||i==1)
                printf("\n");
            else
                printf(" ");
            if((i-1)%10 == 0)
                printf("\n");

            i++;
        }
    }
}
