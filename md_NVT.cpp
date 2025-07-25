#include <iostream>

#include <fstream>

#include <cmath>

#include <random>

#include <vector> // std::vectorを使用



// --- グローバル定数 ---

const double c_pi = 3.141592653589793238462;    // 円周率
const int Ntot = 80;                            // 全粒子数
const int dim = 2;                              // シミュレーションの次元 (2次元)
const int nsteps = 50000;                       // シミュレーションの総ステップ数
const double dt = 0.001;                       // 1ステップあたりの時間刻み幅
const double mass = 1.0;                        // 粒子の質量 (全ての粒子で共通)
const double rho = 1.0;                         // 系の数密度 (粒子数/面積)
const double L[2] = { std::sqrt((double)Ntot/rho), std::sqrt((double)Ntot/rho)}; // 2次元なので面積から一辺を計算
const double T_target = 2.0;  // 目標温度
const double Q = 5.5;        // 熱浴の慣性質量



// --- ポテンシャルパラメータ (型を明記) ---

const double sigmasq[3] = {1.0*1.0, 0.8*0.8, 0.88*0.88};
const double epsilon[3] = {1.0, 1.5, 0.5};



// --- グローバル変数 (状態) ---

static double r[Ntot][dim], v[Ntot][dim], f[Ntot][dim]; // n番の粒子のd次元目(d=0がx軸、d=1がy軸に対応)の位置座標
static double upot[Ntot]; // ポテンシャルエネルギー用配列
static unsigned char spcp[Ntot];  // 各粒子の種類 (species) を示す配列。このコードでは0か1の値
static double LCC[3]; // ポテンシャルカットオフ補正項用
static double zeta = 0.0;     // 熱浴の摩擦係数 ζ

// --- 関数プロトタイプ ---

void make_init_config(void);
void calcforce(void);
void kick(double);
void drift(double);
void wrapall(void);
void update_zeta(double);


int main (void) {

    make_init_config();


    std::cout << "Step,Time,Temperature,KineticEnergy,PotentialEnergy,TotalEnergy" << std::endl;// CSVヘッダーとして出力

    calcforce();

    for (int step = 0; step < nsteps; step++){

        double time = dt*(double)step;

        

         // --- leap flog 法 ---
        kick(0.5*dt);              // 速度を半ステップ進める
        update_zeta(0.5 * dt);     // ζを半ステップ進める

        drift(dt);                 // 位置を1ステップ進める
        wrapall();
        calcforce();               // 新しい位置での力を計算

        update_zeta(0.5 * dt);     // ζを残り半ステップ進める
        kick(0.5*dt);              // 速度を残り半ステップ進める


        if (step % 100 == 0) {          // 100ステップごとに出力


            // 100ステップごとに全ポテンシャルエネルギーを計算して出力する

            double total_upot = 0.0;

            for (int i = 0; i < Ntot; ++i) {

                total_upot += upot[i];

            }

            double total_kin =0.0 ;

            for (int i = 0; i < Ntot; ++i) {

                total_kin += 0.5*mass*(v[i][0]*v[i][0]+v[i][1]*v[i][1]) ;

            }

            const int g = Ntot * dim ;
            double current_temp = 2.0 * total_kin / (double)g; // kB=1と仮定


            std::cout << step << ", " << time << ", " << current_temp << ", " << total_upot + total_kin <<","<< total_kin << "," << total_upot << std::endl;

        }



    }

    return 0;

}



void kick(double delta_t) {

    double invmass = 1.0/mass;

    for (int n=0; n<Ntot; n++) {

        for (int d=0; d<dim; d++) {

            v[n][d] += delta_t * ( f[n][d] * invmass - zeta * v[n][d]);  // dv/dt = F/m - ζdv 

        }

    }

}



void drift(double delta_t) {

    for (int n=0; n<Ntot; n++) {

        for (int d=0; d<dim; d++) {

            r[n][d] += delta_t * v[n][d];  // r(位置) = t(時間) × v(速度)

        }

    }

}

// ζを更新
void update_zeta(double delta_t) {
    // 現在の運動エネルギーを計算
    double current_kin = 0.0;
    for (int n = 0; n < Ntot; ++n) {
        current_kin += 0.5 * mass * (v[n][0] * v[n][0] + v[n][1] * v[n][1]);
    }

    // 目標エネルギーを計算
    const int g = Ntot * dim ; // 系の自由度
    double target_kin = 0.5 * (double)g * T_target; // kB=1と仮定

    // ζを更新 
    zeta += delta_t * ( 2 / Q * (current_kin - target_kin) );  
}

void wrapall(void) {

    for (int n=0; n<Ntot; ++n) {

        for (int d=0; d<dim; d++) {

            while (r[n][d] < 0.0) r[n][d] += L[d];   // r[n][d] < 0.0: 粒子が箱の「下限」を越えた場合。箱の長さL[d]を足すことで、反対側の「上限」から入ってきたとする。

            while (r[n][d] >= L[d]) r[n][d] -= L[d]; // r[n][d] >= L[d]: 粒子が箱の「上限」を越えた場合。箱の長さL[d]を引くことで、反対側の「下限」から入ってきたとする。

        }

    }

}



void calcforce(void) {

    const double rscut = 2.5;

    const double rcsq = rscut*rscut;



    // 力とポテンシャルの初期化

    for (int n=0; n<Ntot; ++n) {

        for(int d=0; d<dim; d++) f[n][d] = 0.0;

        upot[n] = 0.0;

    }



    // カットオフ補正

    for (int ii=0; ii<3; ii++) {    //  U_sihfted(r) = U(r) - U(r_cut) でポテンシャル全体からU(r_cut)だけ引くことで、r > r_cut ではポテンシャルが0、そして r = r_cut でもポテンシャルがちょうど0になり、エネルギーが連続になる。

        double hexsinv = 1.0/(rcsq*rcsq*rcsq);   // rcsqは (r_cut)^2 なので、これは 1/{(r_cut)^2}^3 =1/(r_cut)^6 を計算している。

        LCC[ii] = 4.0*epsilon[ii]*hexsinv*(hexsinv - 1.0);  //  LJポテンシャルの式

    }



    // 全粒子ペアについてループ

    for (int n=0; n<Ntot; ++n) {

        for (int m=n+1; m<Ntot; ++m) { // m=n+1 から始めることで二重計算を避ける

            double dr[dim];  

            double srij = 0.0;

            for (int d=0; d<dim; d++) {

                dr[d] = r[m][d] - r[n][d];  // 2粒子間の距離の差



                // 周期境界条件（最小イメージ規約）// 周期的境界条件というのは平行移動したところも同じとみなせるもので、粒子間距離 dr が最小となるような物を選ぶべきである。

                if (dr[d] >  0.5 * L[d]) dr[d] -= L[d];

                if (dr[d] <= -0.5 * L[d]) dr[d] += L[d];

                srij += dr[d] * dr[d];  // 距離の2乗

            }



            unsigned char spcipair = spcp[m] + spcp[n]; // 粒子nと粒子mの種類（0か1）を足し合わせ

            double ssq = sigmasq[spcipair];   // そのペアに対応する σ^2



            if (srij < rcsq*ssq) {   // 粒子間距離の2乗 srij が、カットオフ距離の2乗 rcsq より小さいとき

                double eps = epsilon[spcipair];   // そのペアに対応する ε

                double srijiv = 1.0/srij;    // 1/r^2

                double sqsovr = ssq*srijiv;  // (σ/r)^2

                double hexsovr = sqsovr*sqsovr*sqsovr;  // (σ/r)^6

                double factor = 24.0*eps*srijiv*hexsovr*(2.0*hexsovr - 1.0); // F(r) = - dU/dr ,factor = F(r) /r

                double pot = 4.0*eps*hexsovr*(hexsovr - 1.0) - LCC[spcipair]; // ぺアのエネルギー



                for (int d=0; d<dim; d++) {

                    f[n][d] -= factor*dr[d];

                    f[m][d] += factor*dr[d]; // 作用反作用の法則

                }

                upot[n] += 0.5*pot;

                upot[m] += 0.5*pot;

            }

        }

    }

}



void make_init_config(void) {

    const double d_min_sq = 0.8*0.8; // 粒子同士がこれ以上近づいてはいけない距離の「2乗」

    int placed = 0;                  // 配置済みの粒子数をカウントする変数,最初は0

   

    // 乱数生成器

    std::random_device rd;

    std::mt19937 gen(rd());

    std::uniform_real_distribution<double> distX(0.0, L[0]);  // 0からボックスのx方向の長さ

    std::uniform_real_distribution<double> distY(0.0, L[1]);  // 0からボックスのy方向の長さ



    while (placed < Ntot) {

        double x_try = distX(gen); // 新しい粒子を配置する候補となる x座標 をランダムに選ぶ

        double y_try = distY(gen); // 新しい粒子を配置する候補となる y座標 をランダムに選ぶ

        bool good = true;



        for (int i = 0; i < placed; ++i) {

            double dx = x_try - r[i][0];   // r[i][0] は、この関数の前のステップですでに配置が成功し、座標が確定したi番目の粒子のx座標を示す。

            double dy = y_try - r[i][1];

            // 周期境界を考慮

            if (dx >  0.5 * L[0]) dx -= L[0];

            if (dx <= -0.5 * L[0]) dx += L[0];

            if (dy >  0.5 * L[1]) dy -= L[1];

            if (dy <= -0.5 * L[1]) dy += L[1];

            double d2 = dx * dx + dy * dy;

            if (d2 < d_min_sq) {

                good = false;

                break;

            }

        }



        if (good) {

            r[placed][0] = x_try;

            r[placed][1] = y_try;


            spcp[placed] = placed % 2; // 粒子種を0と1に割り振る

            placed++;

        }

    }

    if (placed == Ntot) {
    // 全ての粒子が配置されたら、初速度を与える
    std::mt19937 vel_gen(rd());
    std::normal_distribution<double> dist_norm(0.0, std::sqrt(T_target / mass));
    double v_sum[dim] = {0.0, 0.0};
    for (int n = 0; n < Ntot; ++n) {
        for (int d = 0; d < dim; ++d) {
            v[n][d] = dist_norm(gen);
            v_sum[d] += v[n][d];
        }
    }
    for (int d = 0; d < dim; ++d) v_sum[d] /= (double)Ntot;
    for (int n = 0; n < Ntot; ++n) {
        for (int d = 0; d < dim; ++d) {
            v[n][d] -= v_sum[d];
        }
    } 

    }
         
    
}

