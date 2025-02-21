#include<iostream>
using namespace std;
int main(){
    int a,b,ans=1;
    cout<<"Enter a number :";
    cin>>a;
    cout<<"Enter an exponent: ";
    cin>>b;
    for(int i=1;i<=b;i++){
        ans=ans*a;
    }
    cout<<a<<" raised to "<<b<<" is "<<ans<<endl;
}
