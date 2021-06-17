#pragma once
#include <iostream>
#include <math.h>              // mathematical library
#include <cmath> 

namespace Project2 {



	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;
	using namespace std;

	const double ACCURACY = 1.0e-6;

	#ifndef PI 
	#define PI 3.141592653589793238462643
	#endif


	/// <summary>
	/// Summary for MyForm
	/// </summary>
	public ref class MyForm : public System::Windows::Forms::Form
	{
	public:
		MyForm(void)
		{
			InitializeComponent();
			//
			//TODO: Add the constructor code here
			//
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~MyForm()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::TextBox^ r_text;
	private: System::Windows::Forms::TextBox^ q_text;
	private: System::Windows::Forms::TextBox^ sigma_text;
	protected:

	protected:



	private: System::Windows::Forms::TextBox^ S_text;
	private: System::Windows::Forms::TextBox^ X_text;


	protected:

	protected:






	private: System::Windows::Forms::Label^ label1;
	private: System::Windows::Forms::Label^ label2;
	private: System::Windows::Forms::Label^ label3;
	private: System::Windows::Forms::Label^ label4;
	private: System::Windows::Forms::Label^ label5;
	private: System::Windows::Forms::Label^ label6;
	private: System::Windows::Forms::TextBox^ time_text;

	private: System::Windows::Forms::Label^ label7;
	private: System::Windows::Forms::Label^ label8;
	private: System::Windows::Forms::Label^ label9;
	private: System::Windows::Forms::Label^ BAW;
	private: System::Windows::Forms::Label^ BS;
	private: System::Windows::Forms::Button^ count;


	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			this->r_text = (gcnew System::Windows::Forms::TextBox());
			this->q_text = (gcnew System::Windows::Forms::TextBox());
			this->sigma_text = (gcnew System::Windows::Forms::TextBox());
			this->S_text = (gcnew System::Windows::Forms::TextBox());
			this->X_text = (gcnew System::Windows::Forms::TextBox());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->label4 = (gcnew System::Windows::Forms::Label());
			this->label5 = (gcnew System::Windows::Forms::Label());
			this->label6 = (gcnew System::Windows::Forms::Label());
			this->time_text = (gcnew System::Windows::Forms::TextBox());
			this->label7 = (gcnew System::Windows::Forms::Label());
			this->label8 = (gcnew System::Windows::Forms::Label());
			this->label9 = (gcnew System::Windows::Forms::Label());
			this->BAW = (gcnew System::Windows::Forms::Label());
			this->BS = (gcnew System::Windows::Forms::Label());
			this->count = (gcnew System::Windows::Forms::Button());
			this->SuspendLayout();
			// 
			// r_text
			// 
			this->r_text->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 24, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->r_text->Location = System::Drawing::Point(318, 252);
			this->r_text->Name = L"r_text";
			this->r_text->Size = System::Drawing::Size(114, 53);
			this->r_text->TabIndex = 0;
			this->r_text->Text = L"0.08";
			this->r_text->TextAlign = System::Windows::Forms::HorizontalAlignment::Center;
			// 
			// q_text
			// 
			this->q_text->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 24, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->q_text->Location = System::Drawing::Point(693, 96);
			this->q_text->Name = L"q_text";
			this->q_text->Size = System::Drawing::Size(108, 53);
			this->q_text->TabIndex = 1;
			this->q_text->Text = L"-0.04";
			this->q_text->TextAlign = System::Windows::Forms::HorizontalAlignment::Center;
			// 
			// sigma_text
			// 
			this->sigma_text->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 24, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->sigma_text->Location = System::Drawing::Point(693, 178);
			this->sigma_text->Name = L"sigma_text";
			this->sigma_text->Size = System::Drawing::Size(108, 53);
			this->sigma_text->TabIndex = 2;
			this->sigma_text->Text = L"0.2";
			this->sigma_text->TextAlign = System::Windows::Forms::HorizontalAlignment::Center;
			// 
			// S_text
			// 
			this->S_text->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 24, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->S_text->Location = System::Drawing::Point(318, 96);
			this->S_text->Name = L"S_text";
			this->S_text->Size = System::Drawing::Size(114, 53);
			this->S_text->TabIndex = 3;
			this->S_text->Text = L"100";
			this->S_text->TextAlign = System::Windows::Forms::HorizontalAlignment::Center;
			this->S_text->TextChanged += gcnew System::EventHandler(this, &MyForm::S_TextChanged);
			// 
			// X_text
			// 
			this->X_text->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 24, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->X_text->Location = System::Drawing::Point(318, 171);
			this->X_text->Name = L"X_text";
			this->X_text->Size = System::Drawing::Size(114, 53);
			this->X_text->TabIndex = 4;
			this->X_text->Text = L"100";
			this->X_text->TextAlign = System::Windows::Forms::HorizontalAlignment::Center;
			this->X_text->TextChanged += gcnew System::EventHandler(this, &MyForm::X_TextChanged);
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 24, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label1->Location = System::Drawing::Point(22, 103);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(213, 46);
			this->label1->TabIndex = 5;
			this->label1->Text = L"stock price";
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 24, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label2->Location = System::Drawing::Point(22, 178);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(268, 46);
			this->label2->TabIndex = 6;
			this->label2->Text = L"exercise price";
			this->label2->Click += gcnew System::EventHandler(this, &MyForm::label2_Click);
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 24, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label3->Location = System::Drawing::Point(22, 259);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(229, 46);
			this->label3->TabIndex = 7;
			this->label3->Text = L"interest rate";
			this->label3->Click += gcnew System::EventHandler(this, &MyForm::label3_Click);
			// 
			// label4
			// 
			this->label4->AutoSize = true;
			this->label4->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 24, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label4->Location = System::Drawing::Point(457, 103);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(219, 46);
			this->label4->TabIndex = 8;
			this->label4->Text = L"payout rate";
			this->label4->Click += gcnew System::EventHandler(this, &MyForm::label4_Click);
			// 
			// label5
			// 
			this->label5->AutoSize = true;
			this->label5->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 24, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label5->Location = System::Drawing::Point(457, 178);
			this->label5->Name = L"label5";
			this->label5->Size = System::Drawing::Size(162, 46);
			this->label5->TabIndex = 9;
			this->label5->Text = L"volatility";
			// 
			// label6
			// 
			this->label6->AutoSize = true;
			this->label6->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 24, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label6->Location = System::Drawing::Point(457, 259);
			this->label6->Name = L"label6";
			this->label6->Size = System::Drawing::Size(108, 46);
			this->label6->TabIndex = 10;
			this->label6->Text = L"Time";
			// 
			// time_text
			// 
			this->time_text->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 24, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->time_text->Location = System::Drawing::Point(693, 256);
			this->time_text->Name = L"time_text";
			this->time_text->Size = System::Drawing::Size(108, 53);
			this->time_text->TabIndex = 11;
			this->time_text->Text = L"0.25";
			this->time_text->TextAlign = System::Windows::Forms::HorizontalAlignment::Center;
			// 
			// label7
			// 
			this->label7->AutoSize = true;
			this->label7->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 24, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label7->Location = System::Drawing::Point(261, 24);
			this->label7->Name = L"label7";
			this->label7->Size = System::Drawing::Size(389, 46);
			this->label7->TabIndex = 12;
			this->label7->Text = L"American Option Put";
			// 
			// label8
			// 
			this->label8->AutoSize = true;
			this->label8->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 24, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label8->Location = System::Drawing::Point(42, 453);
			this->label8->Name = L"label8";
			this->label8->Size = System::Drawing::Size(694, 46);
			this->label8->TabIndex = 13;
			this->label8->Text = L"Bjerksund Stensland approximation = ";
			// 
			// label9
			// 
			this->label9->AutoSize = true;
			this->label9->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 24, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label9->Location = System::Drawing::Point(22, 392);
			this->label9->Name = L"label9";
			this->label9->Size = System::Drawing::Size(714, 46);
			this->label9->TabIndex = 14;
			this->label9->Text = L"Barone-Adesi Whaley approximation = ";
			// 
			// BAW
			// 
			this->BAW->AutoSize = true;
			this->BAW->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 24, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->BAW->Location = System::Drawing::Point(742, 392);
			this->BAW->Name = L"BAW";
			this->BAW->Size = System::Drawing::Size(42, 46);
			this->BAW->TabIndex = 15;
			this->BAW->Text = L"0";
			this->BAW->Click += gcnew System::EventHandler(this, &MyForm::BAW_Click);
			// 
			// BS
			// 
			this->BS->AutoSize = true;
			this->BS->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 24, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->BS->Location = System::Drawing::Point(742, 453);
			this->BS->Name = L"BS";
			this->BS->Size = System::Drawing::Size(42, 46);
			this->BS->TabIndex = 16;
			this->BS->Text = L"0";
			// 
			// count
			// 
			this->count->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 19.8F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->count->Location = System::Drawing::Point(349, 328);
			this->count->Name = L"count";
			this->count->Size = System::Drawing::Size(197, 50);
			this->count->TabIndex = 17;
			this->count->Text = L"count";
			this->count->UseVisualStyleBackColor = true;
			this->count->Click += gcnew System::EventHandler(this, &MyForm::count_Click);
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(8, 16);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(916, 522);
			this->Controls->Add(this->count);
			this->Controls->Add(this->BS);
			this->Controls->Add(this->BAW);
			this->Controls->Add(this->label9);
			this->Controls->Add(this->label8);
			this->Controls->Add(this->label7);
			this->Controls->Add(this->time_text);
			this->Controls->Add(this->label6);
			this->Controls->Add(this->label5);
			this->Controls->Add(this->label4);
			this->Controls->Add(this->label3);
			this->Controls->Add(this->label2);
			this->Controls->Add(this->label1);
			this->Controls->Add(this->X_text);
			this->Controls->Add(this->S_text);
			this->Controls->Add(this->sigma_text);
			this->Controls->Add(this->q_text);
			this->Controls->Add(this->r_text);
			this->Name = L"MyForm";
			this->Text = L"American Put Calculator";
			this->Load += gcnew System::EventHandler(this, &MyForm::MyForm_Load);
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
	private: System::Void label2_Click(System::Object^ sender, System::EventArgs^ e) {
	}
private: System::Void label3_Click(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void label4_Click(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void BAW_Click(System::Object^ sender, System::EventArgs^ e) {
}


	   double n(const double& z) {  // normal distribution function    
		   return (1.0 / sqrt(2.0 * PI)) * exp(-0.5 * z * z);
	   };

	   double N(const double& z) {
		   if (z > 6.0) { return 1.0; }; // this guards against overflow 
		   if (z < -6.0) { return 0.0; };

		   double b1 = 0.31938153;
		   double b2 = -0.356563782;
		   double b3 = 1.781477937;
		   double b4 = -1.821255978;
		   double b5 = 1.330274429;
		   double p = 0.2316419;
		   double c2 = 0.3989423;

		   double a = fabs(z);
		   double t = 1.0 / (1.0 + a * p);
		   double b = c2 * exp((-z) * (z / 2.0));
		   double n = ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
		   n = 1.0 - b * n;
		   if (z < 0.0) n = 1.0 - n;
		   return n;
	   };

	   double option_price_european_put_payout(const double& S, // spot price
		   const double& K, // Strike (exercise) price,
		   const double& r,  // interest rate
		   const double& q,  // yield on underlying
		   const double& sigma,
		   const double& time) {
		   double sigma_sqr = pow(sigma, 2);
		   double time_sqrt = sqrt(time);
		   double d1 = (log(S / K) + (r - q + 0.5 * sigma_sqr) * time) / (sigma * time_sqrt);
		   double d2 = d1 - (sigma * time_sqrt);
		   double put_price = K * exp(-r * time) * N(-d2) - S * exp(-q * time) * N(-d1);
		   return put_price;
	   };


	   double option_price_american_put_approximated_baw(const double& S,
		   const double& X,
		   const double& r,
		   const double& b,
		   const double& sigma,
		   const double& time) {
		   const double sigma_sqr = sigma * sigma;
		   const double time_sqrt = sqrt(time);
		   const double M = 2.0 * r / sigma_sqr;
		   const double NN = 2.0 * b / sigma_sqr;
		   const double K = 1.0 - exp(-r * time);
		   double q1 = 0.5 * (-(NN - 1) - sqrt(pow((NN - 1), 2.0) + (4.0 * M / K)));

		   // now find initial S value 
		   double q1_inf = 0.5 * (-(NN - 1) - sqrt(pow((NN - 1), 2.0) + 4.0 * M));
		   double S_star_star_inf = X / (1.0 - 1.0 / q1_inf);
		   double h1 = (b * time - 2 * sigma * time_sqrt) * (X / (X - S_star_star_inf));
		   double S_seed = S_star_star_inf + (X - S_star_star_inf) * exp(h1);

		   double Si = S_seed;  // now do Newton iterations to solve for S**
		   int no_iterations = 0;
		   double g = 1;
		   double gprime = 1;
		   while ((fabs(g) > ACCURACY)
			   && (fabs(gprime) > ACCURACY) // to avoid non-convergence
			   && (no_iterations++ < 500)
			   && Si > 0.0) {
			   double p = option_price_european_put_payout(Si, X, r, b, sigma, time);
			   double d1 = (log(Si / X) + (b + 0.5 * sigma_sqr) * time) / (sigma * time_sqrt);
			   g = X - Si - p + (1 - exp((b - r) * time) * N(-d1)) * Si / q1;
			   gprime = (1.0 / q1 - 1.0) * (1.0 - exp((b - r) * time) * N(-d1))
				   + (1.0 / q1) * exp((b - r) * time) * (1.0 / (sigma * time_sqrt)) * n(-d1);
			   Si = Si - (g / gprime);
		   };
		   double S_star_star = Si;
		   if (g > ACCURACY) {
			   S_star_star = S_seed;
		   };  // if not found g**
		   double P = 0;
		   double p = option_price_european_put_payout(S, X, r, b, sigma, time);
		   if (S > S_star_star) {
			   double d1 = (log(S_star_star / X)
				   + (b + 0.5 * sigma_sqr) * time) / (sigma * time_sqrt);
			   double A1 = -(S_star_star / q1) * (1 - exp((b - r) * time) * N(-d1));
			   P = p + A1 * pow((S / S_star_star), q1);
		   }
		   else {
			   P = X - S;
		   };
		   return max(P, p);  // should not be lower than Black Scholes value
	   };

	   inline double phi(double S, double T, double gamma, double H, double X, double r, double b, double sigma) {
		   double sigma_sqr = pow(sigma, 2);
		   double kappa = 2.0 * b / sigma_sqr + 2.0 * gamma - 1.0;
		   double lambda = (-r + gamma * b + 0.5 * gamma * (gamma - 1.0) * sigma_sqr) * T;
		   double d1 = -(log(S / H) + (b + (gamma - 0.5) * sigma_sqr) * T) / (sigma * sqrt(T));
		   double d2 = -(log((X * X) / (S * H)) + (b + (gamma - 0.5) * sigma_sqr) * T) / (sigma * sqrt(T));
		   double phi = exp(lambda) * pow(S, gamma) * (N(d1) - pow((X / S), kappa) * N(d2));
		   return phi;
	   };

	   double option_price_european_call_payout(const double& S, // spot price
		   const double& X, // Strike (exercise) price,
		   const double& r,  // interest rate
		   const double& q,  // yield on underlying
		   const double& sigma, // volatility
		   const double& time) { // time to maturity
		   double sigma_sqr = pow(sigma, 2);
		   double time_sqrt = sqrt(time);
		   double d1 = (log(S / X) + (r - q + 0.5 * sigma_sqr) * time) / (sigma * time_sqrt);
		   double d2 = d1 - (sigma * time_sqrt);
		   double call_price = S * exp(-q * time) * N(d1) - X * exp(-r * time) * N(d2);
		   return call_price;
	   };

	   double option_price_american_call_approximated_bjerksund_stensland(const double& S,
		   const double& K,
		   const double& r,
		   const double& b,
		   const double& sigma,
		   const double& T) {

		   double sigma_sqr = pow(sigma, 2);
		   double B0 = max(K, (r / (r - b) * K));
		   double beta = (0.5 - b / sigma_sqr) + sqrt(pow((b / sigma_sqr - 0.5), 2) + 2.0 * r / sigma_sqr);
		   double Binf = beta / (beta - 1.0) * K;
		   double hT = -(b * T + 2.0 * sigma * sqrt(T)) * ((K * K) / (Binf - B0));
		   double XT = B0 + (Binf - B0) * (1.0 - exp(hT));
		   double alpha = (XT - K) * pow(XT, -beta);
		   double C = alpha * pow(S, beta);
		   C -= alpha * phi(S, T, beta, XT, XT, r, b, sigma);
		   C += phi(S, T, 1, XT, XT, r, b, sigma);
		   C -= phi(S, T, 1, K, XT, r, b, sigma);
		   C -= K * phi(S, T, 0, XT, XT, r, b, sigma);
		   C += K * phi(S, T, 0, K, XT, r, b, sigma);
		   double c = option_price_european_call_payout(S, K, r, b, sigma, T);
		   return max(c, C);


	   };


	   double option_price_american_put_approximated_bjerksund_stensland(const double& S,
		   const double& X,
		   const double& r,
		   const double& q,
		   const double& sigma,
		   const double& T) {

		   return option_price_american_call_approximated_bjerksund_stensland(X, S, r - (r - q), r - q, sigma, T);
	   };


	   double test_baw_approximation_put(const double& S,
		   const double& X,
		   const double& r,
		   const double& q,
		   const double& sigma,
		   const double& T) {

		   //cout << " Put price using Barone-Adesi Whaley approximation = "
			   //<< option_price_american_put_approximated_baw(S, X, r, q, sigma, T) << endl;

		   return option_price_american_put_approximated_baw(S, X, r, q, sigma, T);
	   };

	   void approximations_examples() {
		   //cout << "------------------------------------" << endl;
		   //cout << "Approximations put chapter " << endl;
		   //cout << "------------------------------------" << endl;
		   double S = 100;   double X = 100;     double sigma = 0.20;
		   double r = 0.08;  double b = -0.04;   double time = 0.25;
		   double bs = option_price_american_put_approximated_bjerksund_stensland(S, X, r, b, sigma, time);
		   double baw = test_baw_approximation_put(S, X, r, b, sigma, time);
		   //cout << " Put price using bjerksund_stensland approximation = "
		   //    << bs << endl;

		   
	   };

private: System::Void count_Click(System::Object^ sender, System::EventArgs^ e) {
	
	double S = System::Convert::ToDouble(S_text->Text);   
	double X = System::Convert::ToDouble(X_text->Text);
	double sigma = System::Convert::ToDouble(sigma_text->Text);
	double r = System::Convert::ToDouble(r_text->Text);
	double b = System::Convert::ToDouble(q_text->Text);
	double time = System::Convert::ToDouble(time_text->Text);
	double bs = option_price_american_put_approximated_bjerksund_stensland(S, X, r, b, sigma, time);
	double baw = test_baw_approximation_put(S, X, r, b, sigma, time);

	baw = round(10000 * baw) / 10000;
	bs = round(10000 * bs) / 10000;

	BAW->Text = System::Convert::ToString(baw);
	BS->Text = System::Convert::ToString(bs);

}
private: System::Void S_TextChanged(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void X_TextChanged(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void MyForm_Load(System::Object^ sender, System::EventArgs^ e) {
}
};
}
