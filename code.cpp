#include <bits/stdc++.h>
using namespace std;

struct term {
    int a, b, c, d;
    term(): a(0), b(0), c(0), d(0) {}
    term(int _a, int _b, int _c, int _d): a(_a), b(_b), c(_c), d(_d) {}
    bool operator==(const term& o) const { return b==o.b && c==o.c && d==o.d; }
    bool operator!=(const term& o) const { return !(*this==o); }
};

struct poly {
    int n;
    term *t;
    poly(): n(0), t(NULL) {}
    poly(int _n): n(_n) { t = new term[n]; }
    poly(const poly &p) { n = p.n; t = new term[n]; for (int i=0;i<n;++i) t[i]=p.t[i]; }
    poly& operator=(const poly &o){ if(this==&o) return *this; if(t) delete [] t; n=o.n; t = new term[n]; for(int i=0;i<n;++i) t[i]=o.t[i]; return *this; }
    ~poly(){ if(t) delete [] t; }

    static bool cmp_key(const term &x, const term &y){
        if (x.b != y.b) return x.b > y.b;
        if (x.c != y.c) return x.c > y.c;
        if (x.d != y.d) return x.d > y.d;
        return false;
    }

    void simplify(){
        if (n<=0){ if(t){ delete [] t; t=NULL; } n=1; t=new term[1]; t[0]=term(0,0,0,0); return; }
        sort(t, t+n, cmp_key);
        int m=0; // merged count
        for (int i=0;i<n;){
            term cur = t[i]; long long sum = cur.a;
            int j=i+1; while(j<n && t[j].b==cur.b && t[j].c==cur.c && t[j].d==cur.d){ sum += t[j].a; ++j; }
            if (sum!=0){ t[m++] = term((int)sum, cur.b, cur.c, cur.d); }
            i=j;
        }
        if (m==0){ delete [] t; n=1; t=new term[1]; t[0]=term(0,0,0,0); return; }
        term *nt = new term[m]; for(int i=0;i<m;++i) nt[i]=t[i]; delete [] t; t=nt; n=m;
        // already sorted
    }

    poly operator+(const poly &o) const{
        poly ans(n + o.n);
        for(int i=0;i<n;++i) ans.t[i]=t[i];
        for(int i=0;i<o.n;++i){ ans.t[n+i]=o.t[i]; }
        ans.simplify();
        return ans;
    }
    poly operator-(const poly &o) const{
        poly ans(n + o.n);
        for(int i=0;i<n;++i) ans.t[i]=t[i];
        for(int i=0;i<o.n;++i){ ans.t[n+i]=o.t[i]; ans.t[n+i].a *= -1; }
        ans.simplify();
        return ans;
    }
    poly operator*(const poly &o) const{
        poly ans(n * o.n);
        for(int i=0;i<n;++i){
            for(int j=0;j<o.n;++j){
                long long aa = 1LL*t[i].a * o.t[j].a;
                ans.t[i*o.n + j] = term((int)aa, t[i].b + o.t[j].b, t[i].c + o.t[j].c, t[i].d + o.t[j].d);
            }
        }
        ans.simplify();
        return ans;
    }

    poly derivate() const{
        // derivative of sum of terms
        // for each term a x^b sin^c x cos^d x
        // d/dx = a*b x^{b-1} sin^c cos^d + a*c x^b sin^{c-1} cos^{d+1} - a*d x^b sin^{c+1} cos^{d-1}
        int cnt=0;
        for(int i=0;i<n;++i){ if (t[i].b>0) ++cnt; if (t[i].c>0) ++cnt; if (t[i].d>0) ++cnt; }
        if (cnt==0){ poly z(1); z.t[0]=term(0,0,0,0); return z; }
        poly ans(cnt);
        int k=0;
        for(int i=0;i<n;++i){ const term &u = t[i];
            if (u.b>0){ ans.t[k++] = term(u.a * u.b, u.b-1, u.c, u.d); }
            if (u.c>0){ ans.t[k++] = term(u.a * u.c, u.b, u.c-1, u.d+1); }
            if (u.d>0){ ans.t[k++] = term(-u.a * u.d, u.b, u.c+1, u.d-1); }
        }
        ans.simplify();
        return ans;
    }

    bool is_one() const{ return n==1 && t[0].a==1 && t[0].b==0 && t[0].c==0 && t[0].d==0; }
    bool is_zero() const{ return n==1 && t[0].a==0 && t[0].b==0 && t[0].c==0 && t[0].d==0; }

    string to_string_terms(int &term_count) const{
        term_count = 0;
        for(int i=0;i<n;++i) if (t[i].a!=0) ++term_count;
        if (term_count==0) return string("0");
        string s;
        bool first=true;
        for(int i=0;i<n;++i){ int a = t[i].a; if (a==0) continue; int b=t[i].b, c=t[i].c, d=t[i].d;
            bool isConst = (b==0 && c==0 && d==0);
            if (first){
                if (isConst){ s += std::to_string(a); }
                else{
                    if (a==-1) s.push_back('-');
                    else if (a!=1) s += std::to_string(a);
                    // vars
                    if (b){ s.push_back('x'); if (b>1){ s.push_back('^'); s += std::to_string(b);} }
                    if (c){ s += "sin"; if (c>1){ s.push_back('^'); s += std::to_string(c);} s.push_back('x'); }
                    if (d){ s += "cos"; if (d>1){ s.push_back('^'); s += std::to_string(d);} s.push_back('x'); }
                }
                first=false;
            }else{
                if (a>0) s.push_back('+'); else s.push_back('-');
                int aa = (a>0? a : -a);
                if (isConst){ s += std::to_string(aa); }
                else{
                    if (aa!=1) s += std::to_string(aa);
                    if (b){ s.push_back('x'); if (b>1){ s.push_back('^'); s += std::to_string(b);} }
                    if (c){ s += "sin"; if (c>1){ s.push_back('^'); s += std::to_string(c);} s.push_back('x'); }
                    if (d){ s += "cos"; if (d>1){ s.push_back('^'); s += std::to_string(d);} s.push_back('x'); }
                }
            }
        }
        return s;
    }
};

struct frac{
    poly p, q;
    frac(){}
    frac(int x){ p=poly(1); p.t[0]=term(x,0,0,0); q=poly(1); q.t[0]=term(1,0,0,0);} 
    frac(term _p){ p=poly(1); p.t[0]=_p; q=poly(1); q.t[0]=term(1,0,0,0);} 
    frac(const poly &_p, const poly &_q): p(_p), q(_q) {}

    frac operator+(const frac &o) const{ return frac(p*o.q + q*o.p, q*o.q); }
    frac operator-(const frac &o) const{ return frac(p*o.q - q*o.p, q*o.q); }
    frac operator*(const frac &o) const{ return frac(p*o.p, q*o.q); }
    frac operator/(const frac &o) const{ return frac(p*o.q, q*o.p); }

    frac derivate() const{ return frac(p.derivate()*q - q.derivate()*p, q*q); }

    string to_output_string() const{
        // special cases
        if (p.is_zero()) return string("0");
        if (q.is_one()){ int cnt=0; return p.to_string_terms(cnt); }
        int cntp=0, cntq=0; string sp = p.to_string_terms(cntp); string sq = q.to_string_terms(cntq);
        string out;
        if (cntp>1){ out.push_back('('); out += sp; out.push_back(')'); }
        else{ out += sp; }
        out.push_back('/');
        if (cntq>1){ out.push_back('('); out += sq; out.push_back(')'); }
        else{ out += sq; }
        return out;
    }
};

struct Parser{
    const string &s; int n; int i;
    Parser(const string &str): s(str), n((int)str.size()), i(0) {}

    void skip(){ while(i<n && s[i]==' ') ++i; }

    bool starts_with(const char* lit){ int L=strlen(lit); if (i+L>n) return false; for(int k=0;k<L;++k){ if(s[i+k]!=lit[k]) return false; } return true; }

    int parse_int(){ long long val=0; int start=i; while(i<n && isdigit((unsigned char)s[i])){ val = val*10 + (s[i]-'0'); ++i; }
        if (start==i) return 0; return (int)val; }

    frac parse_expression(){ frac a = parse_term(); skip(); while(i<n && (s[i]=='+' || s[i]=='-')){ char op = s[i++]; frac b = parse_term(); if (op=='+') a = a + b; else a = a - b; skip(); } return a; }

    frac parse_term(){ frac a = parse_factor(); skip(); while(i<n && (s[i]=='*' || s[i]=='/')){ char op = s[i++]; frac b = parse_factor(); if (op=='*') a = a * b; else a = a / b; skip(); } return a; }

    frac parse_factor(){ skip(); if (i<n && s[i]=='+'){ ++i; return parse_factor(); } if (i<n && s[i]=='-'){ ++i; frac f = parse_factor(); return frac(0) - f; } return parse_primary(); }

    frac parse_primary(){ skip(); if (i<n && s[i]=='('){ ++i; frac inside = parse_expression(); skip(); if (i<n && s[i]==')') ++i; return inside; }
        // monomial or integer
        term tm(1,0,0,0); bool saw=false; // coefficient multiplied progressively
        // optional coefficient number (positive)
        if (i<n && isdigit((unsigned char)s[i])){ int num = parse_int(); tm.a = tm.a * num; saw=true; }
        // sequence of tokens x, sin^k x, cos^k x
        while(i<n){ if (s[i]=='x'){
                ++tm.b; ++i; if (i<n && s[i]=='^'){ ++i; int e = parse_int(); tm.b += e-1; }
                saw=true; continue; }
            if (starts_with("sin")){
                i+=3; int e=1; if (i<n && s[i]=='^'){ ++i; e = parse_int(); } if (i<n && s[i]=='x') ++i; tm.c += e; saw=true; continue; }
            if (starts_with("cos")){
                i+=3; int e=1; if (i<n && s[i]=='^'){ ++i; e = parse_int(); } if (i<n && s[i]=='x') ++i; tm.d += e; saw=true; continue; }
            break; }
        if (!saw){ // should not happen on valid input, treat as 0
            return frac(0);
        }
        poly pp(1); pp.t[0]=tm; poly qq(1); qq.t[0]=term(1,0,0,0);
        return frac(pp, qq);
    }
};

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    string str; if(!(cin>>str)) return 0;
    Parser P(str);
    frac g = P.parse_expression();
    // ensure polys are simplified (constructors/ops should already do)
    cout << g.to_output_string() << "\n";
    frac h = g.derivate();
    cout << h.to_output_string() << "\n";
    return 0;
}
