
#include <stdio.h>
#include <stdlib.h>

int Stepen(int base, int exp){
    int i=0, p=1;
    while(++i<exp+1 && (p*=base));
    return p;
}

void MatrixAdd(int* A, int* B, int* C, int n, int x){
    int i,j;
    int m = n/2;
    if (x==0)
    	m = n;
    for (i=0;i<m;i++)
        for (j=0;j<m;j++)
            *(C+i*m+j) = *(A+i*n+j) + *(B+i*n+j);
}

void MatrixSub(int* A, int* B, int* C, int n, int x){
    int i,j;
    int m = n/2;
    if (x==0)
    	m = n;
    for (i=0;i<m;i++)
        for (j=0;j<m;j++)
            *(C+i*m+j) = *(A+i*n+j) - *(B+i*n+j);
}

void Strassen(int* A, int* B, int* C, int n){
    int i,j;
    if(n==2){
        int P=(*A+*(A+n+1))*(*B+*(B+n+1));  //P=(A[0][0]+A[1][1])*(B[0][0]+B[1][1])
        int Q=(*(A+n)+*(A+n+1))*(*B);   //Q=(A[1][0]+A[1][1])*B[0][0]
        int R=(*A)*(*(B+1)-*(B+n+1));   //R=A[0][0]*(B[0][1]-B[1][1])
        int S=(*(A+n+1))*(*(B+n)-*B);   //S=A[1][1]*(B[1][0]-B[0][0])
        int T=(*A+*(A+1))*(*(B+n+1));   //T=(A[0][0]+A[0][1])*B[1][1]
        int U=(*(A+n)-*A)*(*B+*(B+1));  //U=(A[1][0]-A[0][0])*(B[0][0]+B[0][1])
        int V=(*(A+1)-*(A+n+1))*(*(B+n)+*(B+n+1));  //V=(A[0][1]-A[1][1])*(B[1][0]+B[1][1])

        *C=P+S-T+V;
        *(C+1)=R+T;
        *(C+n)=Q+S;
        *(C+n+1)=P+R-Q+U;
    }
    else{
        int m=n/2, x[m][m], y[m][m], o[n][n];
        for(i=0;i<n;i++)
            for(j=0;j<n;j++)
                o[i][j]=0;

        /*P=(A[0][0]+A[1][1])*(B[0][0]+B[1][1])*/
        int P[m][m];
        MatrixAdd(A, A+m*(n+1), x, n, 1);
        MatrixAdd(B, B+m*(n+1), y, n, 1);
        Strassen(x, y, P, m);

        /*Q=(A[1][0]+A[1][1])*B[0][0]*/
        int Q[m][m];
        MatrixAdd(A+m*n, A+m*(n+1), x, n, 1);
        MatrixAdd(B, o, y, n, 1);
        Strassen(x, y, Q, m);

        /*R=A[0][0]*(B[0][1]-B[1][1])*/
        int R[m][m];
        MatrixAdd(A, o, x, n, 1);
        MatrixSub(B+m, B+m*(n+1), y, n, 1);
        Strassen(x, y, R, m);

        /*S=A[1][1]*(B[1][0]-B[0][0])*/
        int S[m][m];
        MatrixAdd(A+m*(n+1), o, x, n, 1);
        MatrixSub(B+m*n, B, y, n, 1);
        Strassen(x, y, S, m);

        /*T=(A[0][0]+A[0][1])*B[1][1]*/
        int T[m][m];
        MatrixAdd(A, A+m, x, n, 1);
        MatrixAdd(B+m*(n+1), o, y, n, 1);
        Strassen(x, y, T, m);

        /*U=(A[1][0]-A[0][0])*(B[0][0]+B[0][1])*/
        int U[m][m];
        MatrixSub(A+m*n, A, x, n, 1);
        MatrixAdd(B, B+m, y, n, 1);
        Strassen(x, y, U, m);

        /*V=(A[0][1]-A[1][1])*(B[1][0]+B[1][1])*/
        int V[m][m];
        MatrixSub(A+m, A+m*(n+1), x, n, 1);
        MatrixAdd(B+m*n, B+m*(n+1), y, n, 1);
        Strassen(x, y, V, m);


        /* 4 части для результата */
        int W[m][m], X[m][m], Y[m][m], Z[m][m];
        MatrixSub(V,T,x,m,0);
        MatrixAdd(S,x,y,m,0);
        MatrixAdd(P,y,W,m,0); // W=P+S-T+V
        MatrixAdd(R,T,X,m,0); // X==R+T
        MatrixAdd(Q,S,Y,m,0); // Y=Q+S
        MatrixSub(U,Q,x,m,0);
        MatrixAdd(R,x,y,m,0);
        MatrixAdd(P,y,Z,m,0); // Z=P+R-Q+U

         /*C = [[P+S-T+V,R+T],
             	   [Q+S,P+R-Q+U]] */
        for (i=0;i<m;i++)
            for (j=0;j<m;j++){
                *(C+i*n+j) = W[i][j]; //C[0][0]=W
                *(C+i*n+j+m) = X[i][j]; //C[0][1]=X
                *(C+(i+m)*n+j) = Y[i][j]; //C[1][0]=Y
                *(C+(i+m)*n+j+m) = Z[i][j]; //C[1][1]=Z
            }
    }
}

void main()
{
    int i,j,k,m,n,n1,n2,n3,n4,o=0;

    /*Ввод первой матрицы*/
    printf("Ввод Количества строк и столбцов первой матрицы:");
    scanf("%d%d",&n1,&n2);
    int A[n1][n2];
    printf("\nВвод первой матрицы:\n");
    for(i=0;i<n1;i++){
        
        for(j=0;j<n2;j++)
            A[i][j] = rand()%3;
    }

    /*Ввод второй матрицы*/
    printf("\nВвод Количества столбцов второй матрицы:");
    scanf("%d",&n3);
    int B[n2][n3];
    printf("\nВвод второй матрицы:\n\n");
    for(i=0;i<n2;i++){
        
        for(j=0;j<n3;j++)
           B[i][j] = rand()%3;
    }

    /* для Алгоритма нужны матрицы размерности степени двойки ,создаем их из исходных добавляя нулевые строки/столбцы*/
    n4 = n1;//находим макс
    if (n1 <= n2)
    	n4 = n2;
    n = n3;
    if(n3 <=n4)
    	n = n4 ;
 
    while(n>(m=Stepen(2,++o)));//ближайшую степень двойки;
    int a[m][m],b[m][m],C[m][m];
    for(i=0;i<m;i++)//создаем квадратные матрицы из нулей
        for(j=0;j<m;j++){
            a[i][j]=0;
            b[i][j]=0;
        }
    for(i=0;i<n1;i++)//добавляем элементы исходных матриц
        for(j=0;j<n2;j++)
            a[i][j]=A[i][j];
    for(i=0;i<n2;i++)
        for(j=0;j<n3;j++)
            b[i][j]=B[i][j];

    /*Первая*/
    printf("\nПервая Матрица:\n");
    for(i=0;i<n1;i++){
        for(j=0;j<n2;j++)
            printf("%d ",a[i][j]);
        printf("\n");
    }
    /*Вторая матрица*/
    printf("\nВторая матрица:\n");
    for(i=0;i<n2;i++){    
        for(j=0;j<n3;j++)
            printf("%d ",b[i][j]);
        printf("\n");
    }

    Strassen(a,b,C,m);    //Вызов умножения по Штрассену .

    /*Результат*/
    printf("\nРеузльтат:\n");
    for(i=0;i<n1;i++){
        for(j=0;j<n3;j++)
            printf("%d ",C[i][j]);
        printf("\n");
    }
}
