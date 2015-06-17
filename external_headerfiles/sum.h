double mean(double A[], double B[], int sizeofa, int sizeofb, int option){

double mean = 0.0;
//option 1: mean(A)
//option 2: mean(B)
//option 3: mean(A*B)
//option 4: mean(A^2)
//option 5: mean(B^2)
//option 6: mean((A-B)^2)

if(option == 1){
for (int i = 0; i<sizeofa; i++){mean = mean + A[i];}
mean = mean/sizeofa;
}

else if(option == 2){
for (int i = 0; i<sizeofb; i++){mean = mean + B[i];}
mean = mean/sizeofb;
}

else if(option == 3 && sizeofa == sizeofb){
for (int i = 0; i<sizeofa; i++){mean = mean + A[i]*B[i];}
mean = mean/sizeofa;
}

else if(option == 4){
for (int i = 0; i<sizeofa; i++){mean = mean + A[i]*A[i];}
mean = mean/sizeofa;
}

else if(option == 5){
for (int i = 0; i<sizeofb; i++){mean = mean + B[i]*B[i];}
mean = mean/sizeofb;
}

else if(option == 6){
for (int i = 0; i<sizeofb; i++){mean = mean + (A[i]-B[i])*(A[i]-B[i]);}

mean = mean/sizeofb;
}
return mean;
}
