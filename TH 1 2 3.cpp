//Nguyen Pham Thanh Hung - 6151071056
//Le Huynh Anh Quy - 6151071090
//Nguyen Huu Tri - 6151071106
#include<stdio.h>
#include<math.h>
#include<conio.h>
#include<cstdlib>
#define MAX 10
void Lagrange(){
    int n , i , j;
    float c , x[100], y[100] , P ,L;
    printf ("Nhap n=");
    scanf ("%d", &n);
    printf("Nhap C=");
    scanf("%f",&c);
    for (i=0 ; i<n ; i++){
        printf ("Nhap x[%d], y[%d]" , i , i );
        scanf ("%f%f", &x[i], &y[i]);
    }
    P=0;
    for (i=0 ; i<n ; i++){
        L=1;
        for(j=0 ; j<n ; j++){
            if(i != j){
                L =L*(c-x[j])/(x[i]-x[j]);
            }
            P=P+y[i]*L;
            printf ("L=%f \n" ,L);

        }
    printf("P(%f)=%f" , c ,P);
    getch();
    }
}
void Nhap(float x[], float y[], int n) {
    for(int i=0; i<n; i++) {
        printf("Nhap moc noi suy thu %d: \n", i);
        printf("Nhap x[%d]= ", i);
        scanf("%f", &x[i]);
        printf("Nhap y[%d]= ", i);
        scanf("%f", &y[i]);
    }
}

void Xuat(float x[], float y[], int n) {
    printf("x: ");
    for(int i=0; i<n; i++) {
        printf("%-10.2f\t", x[i]);
    }
    printf("\n");
    printf("y: ");
    for(int i=0; i<n; i++) {
        printf("%-10.4f\t", y[i]);
    }
    printf("\n");
}

void TinhHieuTy(float y[], float p[100][100], int n) {
    for(int i=0; i<n; i++) {
        p[i][0] = y[i];
    }
    int k;
    for(int j=0; j<n; j++) {
        k = j+1;
        for(int t = 1; t<n; t++) {
            p[t][k] = (p[t][k-1] - p[t-1][k-1]);
        }
    }
}

void BangHieuTy(float p[100][100], int n) {
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            printf("%-10.3f", p[i][j]);
        }
        printf("\n");
    }
}

long tinhGiaithua(int n) {
   if (n > 0) {
       return n * tinhGiaithua(n - 1);
   } 
   else {
        return 1;
   }

}

float Newton(float p[100][100], int n, float c) {
    float kq = 0, t, k = 0;
    int v = 1;
    t = c - 1;
    for(int j = 0; j<n; j++) {
        v = tinhGiaithua(j);
        k += ((t-j)*p[j+1][j+1])/v;
    }
    kq = p[0][0] + k;
    return kq;
}
void NewtonTrenLuoiDeu(){
    int n;
    float x[100], y[100], p[100][100], c, nkd;
    printf("Nhap so moc noi suy: ");
    scanf("%d", &n);
    Nhap(x, y, n);
    printf("Bang noi suy Newton: \n");
    Xuat(x, y, n);
    TinhHieuTy(y, p, n);
    printf("Bang hieu ty: \n");
    BangHieuTy(p, n);
    printf("Nhap gia tri can tinh gan dung: ");
    scanf("%f", &c);
    nkd = Newton(p, n , c);
    printf("Ket qua la: %.5f", nkd);
}
void NewtonTrenLuoiKhongDeu(){
    int n,i,l,j,t,k;
	float x[100], y[100], p[100][100], c, kq;
	//nhap du lieu
	printf(" nhap vao n:");
	scanf("%d",&n);
	for(i=0; i<n;i++){
		printf("nhap vao x[%d], y[%d]", i,i);
		scanf("%f%f", &x[i],&y[i]);
	}
	//tim bang ty hieu
	//gan du lieu x y
	for(i=0; i<n; i++){
		p[i][0] = x[i];
		p[i][1 ]=y[i];
	}
	//tinh ty hieu
	l=1;
	for(j=0; j<n; j++){
		k=j+2;
		for(t=1; t<n; t++){
			p[t][k]=(p[t][k-1]-p[t-1][k-1])/(x[t]-x[t-1]);
		}
		l=l+1;
	}
	//xuat bang ty hieu
	for(i=0; i<n ; i++){
		for(j=0; j<n+2; j++){
			printf("%.2f\t",p[i][j]);
		}
		printf("\n");
	}
	//nhao vao gia tri muon tinh gan dung
	printf("nhap vao gia tri muon tinh gan dung:");
	scanf("%f",&c);
	float tich;
	for(i=0; i<n; i++){
		tich=1;
		for(j=0; j<1; j++){
			tich=tich * (c-x[j]);
		}
		kq=kq+p[i][i+1]*tich;
	}
	printf("ket qua: %.2f", kq);
	getch();
}
enum check { F = 0 , T=1 };
double ff(float x)
{
	double tri ;
	double y;

	y = -x*x ;
	tri = exp(y) ;
	return tri ;
}
void Simpson(){
    int i , n ;
	double a , b , h ,epslon ;
	double s0 , s1 , s2 , i1 , i2;
	enum check t ;
	printf("TICH PHAN DUNG PHUONG PHAP SIMPSON CHO HAM   y = exp(-x*x)\n");
	printf("Can duoi a=");
	fflush(stdin);
	scanf("%lf",&a);
	printf("Can tren b=");
	fflush(stdin);
	scanf("%lf",&b);
	printf("Sai so  epslon=");
	fflush(stdin);
	scanf("%lf",&epslon);
	s2 = 0;
	n = 2 ;
	h = (b - a)/2.0 ;
	s1 = ff(a + h) ;
	s0 = ff(a) + ff(b) ;
	i2 = h*(s0 + 4*s1 + 2*s2)/3.0;
	printf("TICH PHAN SIMPSON CHO HAM   y = exp(-x*x)\n");
	printf("CHO HAM   y = exp(-x*x)\n");
	printf("TREN DOAN [a,b]=[%lf,%lf]\n" , a , b);
	printf("DO LECH epslon = %12.11lf\n" , epslon);
	t = F ;
	do
	{
		i1 = i2 ;
		s2 +=s1 ;
		h /=2 ;
		s1 = 0 ;
		for (i = 1 ; i <= n ;i++) 	s1 +=ff( a + (2*i -1)*h ) ;
		n *= 2 ;
		i2 = h*(s0 + 4*s1 + 2*s2)/3.0 ;
		if ( fabs(i2 - i1) < epslon )
		{
			t = T ;
			printf("\tTICH PHAN    I = %12.3lf\n" ,i2);
		}
	} while (!t) ;
	getch();
}
float FF(float x){
    return x*x;
}
void HinhThang(){
    float a,b;
    int n; 
    printf("\nMoi ban nhap n: ");
    scanf("%d",&n);
    printf("Moi ban nhap a: ");
    scanf("%f",&a);
    printf("\nMoi ban nhap b: ");
    scanf("%f", &b);
    float h=fabs(b-a)/n, x, sum=0;
    int i;
    for(i=1;i<n;i++){
        x=a+i*h;
        sum=sum+FF(x);
    }
    float tichphan = (h/2)*(FF(a)+FF(b)+2*sum);
    printf("\nTich phan gan dung theo phuong phap hinh thang = %f",tichphan);
}

double A[MAX][MAX], B[MAX], X[MAX];
void gauss_dx(int n)
{
  int i=0, j, done=0, m, k;
  double max, c;
  printf("\nTinh nghiem cua he phuong trinh");
  while (!done)
  {
    if (A[i][i] == 0)
    {
      done = 1;
      printf("\nCo phan tu tren duong cheo chinh bang 0");
    }
    if (A[i][i] != 0)
    {
      c = 1/A[i][i];
      for (j=i; j<n; j++)
        A[i][j] = A[i][j] * c;
      B[i] = B[i] * c;
      for (k=i+1; k<n; k++)
      {
        for (j=i+1; j<n; j++)
          A[k][j] = A[k][j] - A[i][j]*A[k][i];
        B[k] = B[k] - B[i] * A[k][i];
        A[k][i] = 0;
      }
    }
    printf("\nLan khu hang %d", i);
    for (k=0; k<n; k++)
    {
      printf("\n");
      for (j=0; j<n; j++)
        printf("%10.5lf", A[k][j]);
      printf(" = %10.5lf", B[k]);
    }
    i++;
    if (i>=n)
      done = 1;
  }
  if (i >= n)
  {
   X[n-1] = B[n-1]/A[n-1][n-1];
   for (i=n-2; i>=0; i--)
   {
     X[i] = 0;
     for (j=n-1; j>i; j--)
       B[i] = B[i] - A[i][j] * X[j];
     X[i] = B[i]/A[i][j];
   }
  }
}

void in_A(int n)
{
  int i, j;
  printf("\nMa tran A :");
  for (i=0; i<n; i++)
  {
    printf("\n");
    for (j=0; j<n; j++)
      printf("%10.5lf", A[i][j]);
  }
}

void in_B(int n)
{
  int i;
  printf("\nMa tran B :\n");
  for (i=0; i<n; i++)
    printf("%10.5lf", B[i]);
}

void in_X(int n)
{
  int i;
  printf("\nMa tran nghiem X :\n");
  for (i=0; i<n; i++)
    printf("%10.5lf", X[i]);
}
void Gauss(){
    int n, i, j;

  printf("Giai he phuong trinh tuyen tinh AX = B.");
  printf("\nbang phuong phap khu GAUSS.");
  printf("\nCho biet cap ma tran : ");
  scanf("%d%*c", &n);
  printf("\nNhap ma tran A :\n");
  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
    {
      printf("A[%d][%d] = ", i, j);
      scanf("%lf", &A[i][j]);
    }
  }
  printf("\nNhap ma tran B :\n");
  for (i=0; i<n; i++)
  {
    printf("B[%d] = ", i);
    scanf("%lf", &B[i]);
  }
  in_A(n);
  in_B(n);
  gauss_dx(n);

  in_X(n);
    
}

void Jacobi(){
    int i , j , k , n , max ;
	float  c , s , epslon;
	enum check t ;
	float b[50] , x0[50] , x1[50] ;
	float a[50][50];

	system("cls");
	printf("PHUONG PHAP LAP JACOBI GIAI HE  \n");
	printf("SO PHUONG TRINH n = ");
	scanf("%d",&n);
	printf("DOC SO LIEU MA TRAN A \n");
	for ( i = 0 ; i < n ; i++ )
	{
		printf("HANG THU %d :" , i);
		for ( j = 0 ; j < n ; j++ )		
			scanf("%f" , &a[i][j]);
	}
	printf("\tMA TRAN A\n");
	for (i = 0 ; i < n ; i++)
	{
		for (j = 0 ; j < n ; j++)	
			printf("%5.2f" , a[i][j]);
			printf("\n");
	}
	printf("DOC VE PHAI B\n");
	for (i = 0 ; i < n ; i++) 
		scanf("%f" , &b[i]);
	printf("\tIN VE PHAI B\n");
	for (i = 0 ; i < n ; i++) 
	{
		printf("%5.2f" , b[i]);
		printf("\n");
	}
	printf("\n");
	printf("DOC VECTOR BAN DAU x0\n");
	for (i = 0 ; i < n ; i++) 
		scanf("%f" , &x0[i]);
	printf("VECTOR BAN DAU x0\n");
	for (i = 0 ; i < n ; i++) 
		printf("%5.2f" , x0[i]);
	printf("\nLAP MA TRAN A*\n");
	for (i = 0 ; i < n ; i++)
	{
		c = 1.0 / a[i][i] ;
		for (j = 0 ; j < n ; j++)
			if (j != i) a[i][j] *=c ;
		b[i] *= c ;
		a[i][i] = 0 ;
	}
	printf("MA TRAN A*\n");
	for (i = 0 ; i < n ; i++)
		for (j = 0 ; j < n ; j++) 	
		printf("%5.2f\n" , a[i][i] );
	printf("SO BUOC LAP NHIEU NHAT max = ");
	scanf("%d" , &max);
	printf("SAI SO DO BUOC CHON TRUOC epslon = ");
	scanf("%f" , &epslon);
	printf("\n");
	printf("SO BUOC LAP CHON max = %d\n" , max);
	printf("SAI SO CHON epslon = %f\n" , epslon);
	printf("THUC HIEN LAP\n");
	k = 1 ;
	t = F ;
	do
	{
		for (i = 0 ; i < n ; i++)
		{
			x1[i] = b[i] ;
			for (j = 0 ; j < n ; j++)
				x1[i] -= a[i][j] * x0[j] ;
		}
		printf("\n");
		printf("SAU LAN LAP THU %d\n" , k);
		for (i = 0 ; i < n ; i++) 	
			printf("%10.4f\n" , x1[i]) ;
		s = 0.0 ;
		for (i = 0 ; i < n ; i++)  
			s += fabs(x1[i] - x0[i]) ;
		if (s >= epslon)
			for (i = 0 ; i < n ; i++) 	
				x0[i] = x1[i] ;
		if (s < epslon)
		{
			t = T ;
			printf("\n");
			printf("PHEP LAP HOI TU SAU %d BUOC\n" , k);
			printf("\tNGHIEM CUA HE\n");
			for ( i = 0 ; i < n ; i++ )
				printf("x[%d] = %f\n" , i , x1[i] );
		}
		k++ ;
		if ( k > max )
		{
			t = T ;
			printf("%12.6f" , s );
			printf("LAP KHONG HOI TU SAU %d \nBUOC\n" , k - 1);
		}
	} 
	while(!t);
	getch();
}
#define epsi 0.00001
void Nhapdathuc(float P[],int &n)
{
	printf("Nhap bac da thuc n=");
	scanf("%d",&n);
	for (int i=n;i>=0;i--)
	{
		printf("Nhap a*x^%.1d=",i);
		scanf("%f",&P[i]);
	}
}
float Tinh(float P[],int n,float x)
{
	float gt=0;
	for (int i=n;i>=0;i--)
	{
		gt=gt+P[i]*pow(x,i);		
	}
	return gt;
}
void ChiaDoi(){
    int n;
	float a[20];
	float x0,x1,x2,y0,y1,y2;
	int maxlap,demlap;
	Nhapdathuc(a,n);
	printf("Tim nghiem phuong trinh bang phuong phap chia doi\n");
	printf("\nCho gia tri x0,x1,maxlap\n");
	printf("Nhap gia tri x0=");
	scanf("%f",&x0);
	printf("\nNhap gia tri x1=");
	scanf("%f",&x1);
	printf("\nNhap so lan max lap: ");
	scanf("%d",&maxlap);
	y0 = Tinh(a,n,x0);
	y1 = Tinh(a,n,x1);
	if((y0*y1)>0)
	{
		printf("\nNghiem khong nam trong khoang x0,x1");
		printf("x0=%.2f\n",x0);
		printf("x1=%.2f\n",x1);
		printf("f(x0)=%.2f\n",y0);
		printf("f(x1)=%.2f\n",y1);
	}	
	demlap=0;
	do
	{
		x2=(x0+x1)/2;
		y2=Tinh(a,n,x2);
		y0=Tinh(a,n,x0);
		if(y0*y2>0)
			x0=x2;
		else
			x1=x2;
		demlap++;
	}
	while (((abs((y2-y0))>epsi)||(demlap<maxlap)));
	if(demlap>maxlap)
	{
		printf("\nPhep lap khong hoi tu sau %d lan lap",maxlap);
		printf("x0=%.2f\n",x0);
		printf("x1=%.2f\n",x1);
		printf("f(x2)=%.2f\n",x2);
	}
	else
	{
		printf("\n\nPhep lap hoi tu sau %d lan lap\n",demlap);
		printf("Nghiem x=%.2f",x2);
	}
	getch();
}

void DayCung(){
    int n;
	float a,b,fa,fb,dx,x;
	float c[20];
	Nhapdathuc(c,n);
	printf("Tim nghiem phuong trinh bang phuong phap day cung\n");
	printf("\nCho gia tri a,b\n");
	printf("Nhap gia tri a=");
	scanf("%f",&a);
	printf("\nNhap gia tri b=");
	scanf("%f",&b);
	fa=Tinh(c,n,a);
	fb=Tinh(c,n,b);
	dx=fa*(b-a)/(fa-fb);
	while(fabs(dx)>epsi)
	{
		x=a+dx;
		fa=Tinh(c,n,x);
		if((fa*fb)<=0)
			a=x;
		else
			b=x;
		fa=Tinh(c,n,a);
		fb=Tinh(c,n,b);
		dx=fa*(b-a)/(fa-fb);
	}
	printf("Nghiem x=%.3f",x);  
    
}
float tinhham (float x)
{
    return x * x * x + x - 5;
}
 
float tinhdaoham(float x)
{
    
    return 3 * x * x + 1;
}
 
void TiepTuyen(){
    float x, y = 0;
    printf("\nNhap x = ");
    scanf("%f", &x);
     
    do {
        y = x;
        float t = tinhham(y) / tinhdaoham(y);
        x = y - t;
        printf("\n%.3f\t%.3f", x, t);
    }
    while (fabs(x - y) > epsi);
     
    printf("\n\nNghiem cua pt la: %f", y);
}
void Lap(){
    float x0,x,q,saiso;
    float f(float);
    printf("Cho sai so = ");
    scanf("%f", &saiso);
    printf("Cho gia tri ban dau cua nghiem = ");
    scanf("%f", &x0);
    x=x0;
    q=f(x);
    if(abs(q-x)>saiso){
        x=q;
        q=f(x);
    }
    printf("Nghiem cua phuong trinh la: %f",q);
    getch();
}
float f(float x){
    float a=exp(1/3)*log(1000-x);
    return a;
}
// lap don
void LapDon(){
    int n;
    float x,x0;
    float f(float);
    printf("Cho so lan lap n= ");
    scanf("%d",&n);
    printf("Cho gia tri ban dau cua nghiem x0= ");
    scanf("%f",&x0);
    x=x0;
    for(int i=1;i<=n;i++)x=f(x);
    printf("Nghiem cua phuong trinh la:%.4f",x);
    getch();
    
}
int main(){
    int n;
    do{
        
           printf("\n\t\t=================================================================\n");
           printf("\t\t=\t\t\t________\t\t\t\t=\n");
           printf("\t\t=\t\t\t  MENU\t\t\t\t\t=\n");
           printf("\t\t=\t\t\t________\t\t\t\t=\n");
           printf("\t\t=\t\t\t\t\t\t\t\t=\n");
           printf("\t\t=\t\t\t1 - TH1\t\t\t\t\t=\n");
           printf("\t\t=\t\t\t2 - TH2\t\t\t\t\t=\n");
           printf("\t\t=\t\t\t3 - TH3\t\t\t\t\t=\n");
           printf("\t\t=\t\t\tx- Thoat\t\t\t\t=\n");
           printf("\t\t=================================================================\n");
           printf("Moi ban nhap lua chon: \n");
           scanf(" %d", &n);
           if(n==1){
               char gt;
               do{
        
                   printf("\n\t\t=================================================================\n");
                   printf("\t\t=\t\t\t________\t\t\t\t=\n");
                   printf("\t\t=\t\t\t MENU 1\t\t\t\t\t=\n");
                   printf("\t\t=\t\t\t________\t\t\t\t=\n");
                   printf("\t\t=\t\t\t\t\t\t\t\t=\n");
                   printf("\t\t=\t\t\t1 - Lagrange\t\t\t\t=\n");
                   printf("\t\t=\t\t\t2 - Newton tren luoi deu\t\t=\n");
                   printf("\t\t=\t\t\t3 - Newton tren luoi khong deu\t\t=\n");
                   printf("\t\t=================================================================\n");
                   printf("Moi ban nhap lua chon: \n");
                   scanf(" %c", &gt);
                   if(gt=='1'){
                        Lagrange();
                   }else if(gt=='2'){
                        NewtonTrenLuoiDeu();
                   }else if(gt=='3'){
                        NewtonTrenLuoiKhongDeu();
                   }
                }while(gt!='x');
            }
            if(n==2){
               char gt;
               do{
                   printf("\n\t\t=================================================================\n");
                   printf("\t\t=\t\t\t________\t\t\t\t=\n");
                   printf("\t\t=\t\t\t MENU 2\t\t\t\t\t=\n");
                   printf("\t\t=\t\t\t________\t\t\t\t=\n");
                   printf("\t\t=\t\t\t\t\t\t\t\t=\n");
                   printf("\t\t=\t\t\t1 - Simpson\t\t\t\t=\n");
                   printf("\t\t=\t\t\t2 - Hinh thang\t\t\t\t=\n");
                   printf("\t\t=\t\t\t3 - Gauss\t\t\t\t=\n");
                   printf("\t\t=\t\t\t4 - Jacobi\t\t\t\t=\n");
                   printf("\t\t=================================================================\n");
                   printf("Moi ban nhap lua chon: \n");
                   scanf(" %c", &gt);
                   if(gt=='1'){
                        Simpson();
                   }else if(gt=='2'){
                        HinhThang();
                   }else if(gt=='3'){
                        Gauss();
                   }else if(gt=='4'){
                        Jacobi();
                   }
                }while(gt!='x');
            }
            else if(n==3){
               char gt;
                do{
        
                   printf("\n\t\t=================================================================\n");
                   printf("\t\t=\t\t\t________\t\t\t\t=\n");
                   printf("\t\t=\t\t\t MENU 3\t\t\t\t\t=\n");
                   printf("\t\t=\t\t\t________\t\t\t\t=\n");
                   printf("\t\t=\t\t\t\t\t\t\t\t=\n");
                   printf("\t\t=\t\t\t1 - Chia doi\t\t\t\t=\n");
                   printf("\t\t=\t\t\t2 - Day cung\t\t\t\t=\n");
                   printf("\t\t=\t\t\t3 - Tiep tuyen\t\t\t\t=\n");
                   printf("\t\t=\t\t\t4 - Lap\t\t\t\t\t=\n");
                   printf("\t\t=\t\t\t5 - LapDon\t\t\t\t=\n");
                   printf("\t\t=================================================================\n");
                   printf("Moi ban nhap lua chon: \n");
                   scanf(" %c", &gt);
                    if(gt=='1'){
                        ChiaDoi();
                   }else if(gt=='2'){
                        DayCung();
                   }else if(gt=='3'){
                        TiepTuyen();
                   }else if(gt=='4'){
                        Lap();
                   }else if(gt=='5'){
                        LapDon();
                   }
                }while(gt!='x');
            }
    }while(n>=1 && n<=3);
}
