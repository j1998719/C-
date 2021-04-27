#pragma once

#include <cmath>                 // standard C mathematical library
#include <algorithm>             // defines the max() operator
#include <vector>                // STL vector templates
#include <iostream>


namespace Project2 {

	
	using namespace std;
	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

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
	private: System::Windows::Forms::Label^ welcome_label;
	private: System::Windows::Forms::Button^ count;
	private: System::Windows::Forms::TextBox^ stock_price;
	private: System::Windows::Forms::TextBox^ exec_rate;

	protected:




	private: System::Windows::Forms::Label^ label1;
	private: System::Windows::Forms::Label^ label2;
	private: System::Windows::Forms::Label^ label3;
	private: System::Windows::Forms::TextBox^ r_rate;

	private: System::Windows::Forms::Label^ label4;
	private: System::Windows::Forms::TextBox^ Steps;


	private: System::Windows::Forms::Label^ label5;
	private: System::Windows::Forms::TextBox^ Sigma;

	private: System::Windows::Forms::Label^ label6;
	private: System::Windows::Forms::TextBox^ Time;

	private: System::Windows::Forms::Label^ label7;
	private: System::Windows::Forms::TextBox^ q_rate;
	private: System::Windows::Forms::Label^ result_text;
	private: System::Windows::Forms::Label^ label8;

















	protected:

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
			this->welcome_label = (gcnew System::Windows::Forms::Label());
			this->count = (gcnew System::Windows::Forms::Button());
			this->stock_price = (gcnew System::Windows::Forms::TextBox());
			this->exec_rate = (gcnew System::Windows::Forms::TextBox());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->r_rate = (gcnew System::Windows::Forms::TextBox());
			this->label4 = (gcnew System::Windows::Forms::Label());
			this->Steps = (gcnew System::Windows::Forms::TextBox());
			this->label5 = (gcnew System::Windows::Forms::Label());
			this->Sigma = (gcnew System::Windows::Forms::TextBox());
			this->label6 = (gcnew System::Windows::Forms::Label());
			this->Time = (gcnew System::Windows::Forms::TextBox());
			this->label7 = (gcnew System::Windows::Forms::Label());
			this->q_rate = (gcnew System::Windows::Forms::TextBox());
			this->result_text = (gcnew System::Windows::Forms::Label());
			this->label8 = (gcnew System::Windows::Forms::Label());
			this->SuspendLayout();
			// 
			// welcome_label
			// 
			this->welcome_label->AutoSize = true;
			this->welcome_label->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 19.8F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->welcome_label->Location = System::Drawing::Point(202, 31);
			this->welcome_label->Name = L"welcome_label";
			this->welcome_label->Size = System::Drawing::Size(381, 38);
			this->welcome_label->TabIndex = 0;
			this->welcome_label->Text = L"Bermudan Put Calculator";
			this->welcome_label->Click += gcnew System::EventHandler(this, &MyForm::label1_Click);
			// 
			// count
			// 
			this->count->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 19.8F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->count->Location = System::Drawing::Point(605, 175);
			this->count->Name = L"count";
			this->count->Size = System::Drawing::Size(154, 54);
			this->count->TabIndex = 1;
			this->count->Text = L"Count";
			this->count->UseVisualStyleBackColor = true;
			this->count->Click += gcnew System::EventHandler(this, &MyForm::label1_Click);
			// 
			// stock_price
			// 
			this->stock_price->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 19.8F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->stock_price->Location = System::Drawing::Point(199, 103);
			this->stock_price->Name = L"stock_price";
			this->stock_price->Size = System::Drawing::Size(120, 45);
			this->stock_price->TabIndex = 2;
			this->stock_price->Text = L"80";
			this->stock_price->TextChanged += gcnew System::EventHandler(this, &MyForm::textBox1_TextChanged);
			// 
			// exec_rate
			// 
			this->exec_rate->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 19.8F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->exec_rate->Location = System::Drawing::Point(199, 175);
			this->exec_rate->Name = L"exec_rate";
			this->exec_rate->Size = System::Drawing::Size(122, 45);
			this->exec_rate->TabIndex = 3;
			this->exec_rate->Text = L"100";
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 19.8F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label1->Location = System::Drawing::Point(121, 243);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(28, 38);
			this->label1->TabIndex = 4;
			this->label1->Text = L"r";
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 19.8F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label2->Location = System::Drawing::Point(121, 175);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(39, 38);
			this->label2->TabIndex = 5;
			this->label2->Text = L"K";
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 19.8F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label3->Location = System::Drawing::Point(121, 110);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(39, 38);
			this->label3->TabIndex = 6;
			this->label3->Text = L"S";
			// 
			// r_rate
			// 
			this->r_rate->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 19.8F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->r_rate->Location = System::Drawing::Point(199, 243);
			this->r_rate->Name = L"r_rate";
			this->r_rate->Size = System::Drawing::Size(122, 45);
			this->r_rate->TabIndex = 7;
			this->r_rate->Text = L"0.2";
			// 
			// label4
			// 
			this->label4->AutoSize = true;
			this->label4->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 19.8F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label4->Location = System::Drawing::Point(330, 245);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(101, 38);
			this->label4->TabIndex = 9;
			this->label4->Text = L"Steps";
			// 
			// Steps
			// 
			this->Steps->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 19.8F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->Steps->Location = System::Drawing::Point(451, 243);
			this->Steps->Name = L"Steps";
			this->Steps->Size = System::Drawing::Size(120, 45);
			this->Steps->TabIndex = 8;
			this->Steps->Text = L"500";
			// 
			// label5
			// 
			this->label5->AutoSize = true;
			this->label5->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 19.8F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label5->Location = System::Drawing::Point(327, 178);
			this->label5->Name = L"label5";
			this->label5->Size = System::Drawing::Size(109, 38);
			this->label5->TabIndex = 11;
			this->label5->Text = L"Sigma";
			// 
			// Sigma
			// 
			this->Sigma->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 19.8F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->Sigma->Location = System::Drawing::Point(451, 175);
			this->Sigma->Name = L"Sigma";
			this->Sigma->Size = System::Drawing::Size(120, 45);
			this->Sigma->TabIndex = 10;
			this->Sigma->Text = L"0.25";
			this->Sigma->TextChanged += gcnew System::EventHandler(this, &MyForm::Sigma_TextChanged);
			// 
			// label6
			// 
			this->label6->AutoSize = true;
			this->label6->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 19.8F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label6->Location = System::Drawing::Point(335, 106);
			this->label6->Name = L"label6";
			this->label6->Size = System::Drawing::Size(89, 38);
			this->label6->TabIndex = 13;
			this->label6->Text = L"Time";
			this->label6->Click += gcnew System::EventHandler(this, &MyForm::label6_Click);
			// 
			// Time
			// 
			this->Time->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 19.8F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->Time->Location = System::Drawing::Point(451, 99);
			this->Time->Name = L"Time";
			this->Time->Size = System::Drawing::Size(120, 45);
			this->Time->TabIndex = 12;
			this->Time->Text = L"1.0";
			// 
			// label7
			// 
			this->label7->AutoSize = true;
			this->label7->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 19.8F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label7->Location = System::Drawing::Point(121, 319);
			this->label7->Name = L"label7";
			this->label7->Size = System::Drawing::Size(35, 38);
			this->label7->TabIndex = 15;
			this->label7->Text = L"q";
			this->label7->Click += gcnew System::EventHandler(this, &MyForm::label7_Click);
			// 
			// q_rate
			// 
			this->q_rate->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 19.8F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->q_rate->Location = System::Drawing::Point(199, 312);
			this->q_rate->Name = L"q_rate";
			this->q_rate->Size = System::Drawing::Size(120, 45);
			this->q_rate->TabIndex = 14;
			this->q_rate->Text = L"0.0";
			this->q_rate->TextChanged += gcnew System::EventHandler(this, &MyForm::textBox7_TextChanged);
			// 
			// result_text
			// 
			this->result_text->AutoSize = true;
			this->result_text->BackColor = System::Drawing::SystemColors::AppWorkspace;
			this->result_text->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 19.8F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->result_text->Location = System::Drawing::Point(415, 401);
			this->result_text->MinimumSize = System::Drawing::Size(100, 20);
			this->result_text->Name = L"result_text";
			this->result_text->Size = System::Drawing::Size(100, 38);
			this->result_text->TabIndex = 16;
			// 
			// label8
			// 
			this->label8->AutoSize = true;
			this->label8->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 19.8F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(0)));
			this->label8->Location = System::Drawing::Point(222, 401);
			this->label8->Name = L"label8";
			this->label8->Size = System::Drawing::Size(178, 38);
			this->label8->TabIndex = 17;
			this->label8->Text = L"Put Price =";
			this->label8->Click += gcnew System::EventHandler(this, &MyForm::label8_Click);
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(8, 16);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(784, 486);
			this->Controls->Add(this->label8);
			this->Controls->Add(this->result_text);
			this->Controls->Add(this->label7);
			this->Controls->Add(this->q_rate);
			this->Controls->Add(this->label6);
			this->Controls->Add(this->Time);
			this->Controls->Add(this->label5);
			this->Controls->Add(this->Sigma);
			this->Controls->Add(this->label4);
			this->Controls->Add(this->Steps);
			this->Controls->Add(this->r_rate);
			this->Controls->Add(this->label3);
			this->Controls->Add(this->label2);
			this->Controls->Add(this->label1);
			this->Controls->Add(this->exec_rate);
			this->Controls->Add(this->stock_price);
			this->Controls->Add(this->count);
			this->Controls->Add(this->welcome_label);
			this->Name = L"MyForm";
			this->Text = L"MyForm";
			this->Load += gcnew System::EventHandler(this, &MyForm::MyForm_Load);
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
	private: System::Void MyForm_Load(System::Object^ sender, System::EventArgs^ e) {
	}


		   double option_price_put_bermudan_binomial(const double& S,
			   const double& X,
			   const double& r,
			   const double& q,
			   const double& sigma,
			   const double& time,
			   const vector<double>& potential_exercise_times,
			   const int& steps) {
			   double delta_t = time / steps;
			   double R = exp(r * delta_t);
			   double Rinv = 1.0 / R;
			   double u = exp(sigma * sqrt(delta_t));
			   double uu = u * u;
			   double d = 1.0 / u;
			   double p_up = (exp((r - q) * delta_t) - d) / (u - d);
			   double p_down = 1.0 - p_up;
			   vector<double> prices(steps + 1);
			   vector<double> put_values(steps + 1);

			   vector<int> potential_exercise_steps; // create list of steps at which exercise may happen
			   for (int i = 0; i < potential_exercise_times.size(); ++i) {
				   double t = potential_exercise_times[i];
				   if ((t > 0.0) && (t < time)) {
					   potential_exercise_steps.push_back(int(t / delta_t));
				   };
			   };

			   prices[0] = S * pow(d, steps);  // fill in the endnodes.
			   for (int i = 1; i <= steps; ++i) prices[i] = uu * prices[i - 1];
			   for (int i = 0; i <= steps; ++i) put_values[i] = max(0.0, (X - prices[i])); // put payoffs at maturity
			   for (int step = steps - 1; step >= 0; --step) {
				   bool check_exercise_this_step = false;
				   for (int j = 0; j < potential_exercise_steps.size(); ++j) {
					   if (step == potential_exercise_steps[j]) { check_exercise_this_step = true; };
				   };
				   for (int i = 0; i <= step; ++i) {
					   put_values[i] = (p_up * put_values[i + 1] + p_down * put_values[i]) * Rinv;
					   prices[i] = d * prices[i + 1];
					   if (check_exercise_this_step) put_values[i] = max(put_values[i], X - prices[i]);
				   };
			   };
			   return put_values[0];
		   };

		   double test_bermudan_option(const double& S,
			   const double& K,
			   const double& r,
			   const double& time,
			   const double& sigma,
			   const int& steps,
			   const double& q) {

			   vector<double> potential_exercise_times;  potential_exercise_times.push_back(0.25);
			   potential_exercise_times.push_back(0.5);  potential_exercise_times.push_back(0.75);
			   cout << " Bermudan put price = "
				   << option_price_put_bermudan_binomial(S, K, r, q, sigma, time, potential_exercise_times, steps)
				   << endl;

			   return option_price_put_bermudan_binomial(S, K, r, q, sigma, time, potential_exercise_times, steps);
		   };



	private: System::Void label1_Click(System::Object^ sender, System::EventArgs^ e) {
		double S = System::Convert::ToDouble(stock_price->Text);
		double K = System::Convert::ToDouble(exec_rate->Text);
		double r = ::System::Convert::ToDouble(r_rate->Text);
		double time = System::Convert::ToDouble(Time->Text);
		double sigma = System::Convert::ToDouble(Sigma->Text);
		int steps = System::Convert::ToInt16(Steps->Text);;
		double q = System::Convert::ToDouble(q_rate->Text);


		double bermudan_option = test_bermudan_option(S,K,r,time,sigma,steps,q);
		result_text->Text = System::Convert::ToString(bermudan_option);

	}
	private: System::Void textBox1_TextChanged(System::Object^ sender, System::EventArgs^ e) {
	}

private: System::Void label1_Click_1(System::Object^ sender, System::EventArgs^ e) {
}



private: System::Void label6_Click(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void label7_Click(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void textBox7_TextChanged(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void Sigma_TextChanged(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void label8_Click(System::Object^ sender, System::EventArgs^ e) {
}
};
}

