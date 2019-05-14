#include<iostream>
#include<cmath>
#include <iomanip>
using namespace std;

double Dr=1.22, Tm=1023.15, XA0=0.98, WB=6900,FH2O=0.034, T0=898, P=121, pc=1440, dHr=139000, MA=106, MB=104, Cp=2.177, pai=3.1416, dl=pow(10,-6);
double FA0=WB/MB/24/3600/XA0, Ac=pai*pow(Dr,2)/4, dWc=pc*Ac*dl;//����Ħ�����ʣ���������΢Ԫ���ڴ߻���������

double ra(double T, double XA);//��Ӧ���ʼ���

int main()
{
	double l=dl/10, XA=0, T=T0, L0=0.5, L=0;
	int i=0;
	cout.setf(ios::showpoint);
	cout.precision(6);
	cout.setf(ios::fixed);
	cout << '\t'<< "����߶�l/m" << '\t'<< "��Ӧ�¶�T/K"  << '\t' << "ת����XA " << endl;
	cout << '\t'<< l << '\t' << T << '\t'<< XA << endl;

	while (XA<XA0)
	{
		double dXA, dT,dTd;
		double U;
		U=140*pow(l,-0.33)*4.186/60/1000;
		dTd=((-dHr)*(-ra(T,XA))*dWc-pai*Dr*U*(T-Tm)*dl)/(FA0 * 106 + 18 * FH2O) / Cp;//�����¶ȱ仯
		dXA=(-ra(T,XA))*dWc/FA0;//ת���ʱ仯
		double rA;
		rA=(ra(T+dTd, XA+dXA)+ra(T, XA))/2;
		U=140*pow(l,-0.33)*4.186/60/1000;
		dT=((-dHr)*(-rA)*dWc-pai*Dr*U*(T+dTd/2-Tm)*dl)/(FA0 * 106 + 18 * FH2O) / Cp;

		while(fabs(dTd - dT) > 1E-5)//�����õ�dXA��dT
		{
			dTd=dT;
			dXA=(-rA)*dWc/FA0;;
			rA=(ra(T+dTd, XA+dXA)+ra(T, XA))/2;
			U=140*pow(l,-0.33)*4.186/60/1000;
			dT=((-dHr)*(-rA)*dWc-pai*Dr*U*(T+dTd/2-Tm)*dl)/(FA0 * 106 + 18 * FH2O) / Cp;
		}

		l = l + dl;
		XA = XA+ dXA;
		T = T + dT;

		if (l > L0)//ÿ��0.3m���һ�δ����¶ȡ�ת����
		{
			L0= L0 + 0.3;
			cout << '\t'<< l << '\t' << T << '\t' <<XA  << endl;
		}
	}

	L = l;
	cout << '\t'<< l << '\t' << T << '\t' << XA << endl;
	cout<<endl;
	cout << '\t'<< "�����ܸ�L��" << L<< "m"<<endl;

}

double ra(double T, double XA)//��Ӧ���ʼ���
{
	double k1,k2;
	k1=2.76*pow(10,-6)*exp(-10983 / T + 9.44);
	k2=2.76*pow(10,-6)*exp(-3676.394 / T - 10.525);
	double yA,yB;
	yA = (1 - XA)*FA0 / ((1 + XA)*FA0 + FH2O);
	yB = XA *FA0 / ((1 + XA)*FA0 + FH2O);
	double rA;
	rA=k2*pow(P*yB,2)-k1*P*yA;
	return rA;
}
