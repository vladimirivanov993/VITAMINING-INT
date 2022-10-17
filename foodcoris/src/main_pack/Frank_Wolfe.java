package main_pack;

public class Frank_Wolfe {

	public Frank_Wolfe() {
	}

/**умножение числа на вектор*/	
public static double[] mul(double a, double[] m){
	double[] r=new double[m.length];
	for(int i=0;i<m.length;i++){
		r[i]=a*m[i];
	}
	return r;
	
}

/**стандартное скалярное произведение векторов*/
public static double mul(double[] a, double[] b){
	double r=0;
	for(int i=0;i<Math.min(a.length,b.length);i++){
		r+=a[i]*b[i];
	}
	return r;
}

/**умножение матрицы на вектор*/
public static double[] mul(double[][] A, double[] b){
	double[] r=new double[A.length];
	for(int i=0;i<A.length;i++){
		r[i]=mul(A[i], b);
	}
	return r;
}

/**разность векторов*/
public static double[] minus(double[] a, double[] b){
double[] r=new double[Math.min(a.length,b.length)];
for(int i=0;i<Math.min(a.length, b.length);i++){
r[i]=a[i]-b[i];	
}
return r;
}

/**сумма векторов*/
public static double[] plus(double[] a, double[] b){
double[] r=new double[Math.min(a.length,b.length)];
for(int i=0;i<Math.min(a.length, b.length);i++){
r[i]=a[i]+b[i];	
}
return r;
}


/** градиент <Ax,x>+<b,x> = 2A(uk)+b */ 
public static double[] gradient(double[][] A, double[] b, double[] uk){
	double[] r=mul(A,uk);
	return r=plus(mul(2,r),b);
}



/**умножение числа на матрицу*/
public static double[][] mul(double a, double[][] A){
	for(int i=0;i<A.length;i++){
		for(int j=0;j<A[0].length;j++){
			A[i][j]=a*A[i][j];
		}
	}
	return A;
}


public static void toPrint(double[] a){
//System.out.println();
	for(int i=0;i<a.length;i++){
		System.out.print(a[i]+" ");
	}
	System.out.println();	
}


/**функция <Ax,x>+<b,x>*/
public static double quadrf(double[][] A,double[] b,double[] uk){
	return mul(mul(A,uk),uk)+mul(b,uk); 
}

/**получение решения с использованием симплекс-метода*/
public static double[] getSolSimplex(double[] uk, double[][] A,double[] b, double[][] simp_A, double[] simp_b, double eps){
	double[] s=new Simplex(simp_A, simp_b, gradient(A, b, uk)).primal();
	int k=0;
	double gamma=1;
	while(true){
		gamma=2/(k+2);
		uk=plus(mul(1-gamma,uk), mul(gamma,s));k++; 
		s=new Simplex(simp_A, simp_b, gradient(A, b, uk)).primal();
		toPrint(uk);System.out.println("fx="+(mul(mul(A,uk),uk)+mul(b,uk)  ));
		toPrint(s);System.out.println("fx="+(mul(mul(A,s),s)+mul(b,s)  ));
		if(Math.abs(quadrf(A,b,uk)-quadrf(A,b,s))<eps)break;
		
	}
	return s;
}



/**решение задачи методом перебора точек R2 с шагом eps*/
public static double[] gM(double[][] A, double[] b, double[] c, double eps, double[][] A2, double[] b2){
	double max=Double.NEGATIVE_INFINITY;
	double curMax=0;
	double[] r=new double[3];
	for(double i=0;i<c[1];i+=eps){
		for(double j=0;j<c[2];j+=eps){
			curMax=mul(mul(A,new double[] {i,j}),new double[] {i,j})+mul(b,new double[] {i,j});
		if (curMax>=max&&positive(minus(b2,(mul(A2,new double[]{i,j}))))) {max=curMax;r[0]=curMax;r[1]=i;r[2]=j;}
		}
	}
	return r;
}


/**проверка, что вектор находится в положительном ортанте пространства Rn*/
private static boolean positive(double[] minus) {
	boolean b=true;
	for(int i=0;i<minus.length;i++){	
	b=b&minus[i]>0;	
	}
	return b;
}



public static void main(String[] args){

//double[] uk={0.0,0.0};
double[][] A={ {1,-1},{3,1} };
double[] b={1,7};
double[] c={32,8};
Simplex s=new Simplex(A, b, c);
double[] sp=s.primal();
System.out.print("Проверка симплекс-метода(линейные функция и ограничения): ");toPrint(sp);

double[][] gA={{-1,0},{0,-1}};
double[] gb={32,8};
double[]d=gradient(gA,gb,new double[]{0,0});
Simplex s2=new Simplex(A,b,d);
double[] sp2=s2.primal();
System.out.print("Проверка корректности применения симплекс метода к градиенту квадратичной функции: ");toPrint(sp2);

double[][] gA2={ {4,8},{8,10} };
double[] gb2={7,11};
double[] uk2={2,3};
double[][] A2={{-1,1},{0.75, -1},{0.5,1}};
double[] b2={3,1.5,6};

System.out.print("Проверка корректности функции вычисления градиента { {4,8},{8,10} }+{7,11} в точке: {2,3}: ");toPrint(gradient(gA2,gb2,uk2));
System.out.print("Решение задачи методом условного градиента с функцией (1): ");toPrint(getSolSimplex(uk2, gA2, gb2,  A2, b2,0.01));
System.out.print("Решение задачи методом перебора с функцией (1): ");toPrint(gM(gA2,  gb2, new double[]{0,6,6}, 0.1, A2, b2)); //eps

double[][] gA3={{-40,-8},{-8,-110}};
double[] gb3={700,111};

System.out.print("Решение задачи методом условного градиента с функцией (2): ");toPrint(getSolSimplex(uk2, gA3, gb3,  A2, b2,0.01));
System.out.print("Решение задачи методом перебора с функцией (2): ");toPrint(gM(gA3,  gb3, new double[]{0,6,6}, 0.1, A2, b2)); //eps
}	
}