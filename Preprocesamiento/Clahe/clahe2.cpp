#include "dicom/DicomReader.h"
#include <vector>       // std::vector
#include <cmath>

using namespace std;
    vector<vector<int>> newMatriz(DicomReader dicomObj){
        vector<vector<int>> data = dicomObj.getIntImageMatrix(12); //guarda la matriz de la img dicom
        return data;
    }

    vector<vector<int>> newCLAHE_GO(vector<vector<int>> data){
        vector<vector<int>> CLAHE_GO(data.size());
        copy(data.begin(), data.end(), CLAHE_GO.begin());
        return CLAHE_GO;
    }
    
    vector<vector<int>> claheGo(vector<vector<int>> CLAHE_GO,vector<vector<int>> data,DicomReader dicomObj){
        int _step = 8;
        int block = _step;//pblock
        int width = (sizeof(data)/sizeof(data[0]));
        int height= (sizeof(data[0])/sizeof(data[0][0]));
        int width_block = width/block;
        int height_block = height/block;
        int tmp2[8*8][256] ={{0}};
        float C2[8*8][256] = {{0.0}};

        int total = width_block * height_block;

        for (int i=0;i<block;i++)
        {
            for (int j=0;j<block;j++)
            {

                int start_x = i*width_block;
                int end_x = start_x + width_block;
                int start_y = j*height_block;
                int end_y = start_y + height_block;
                int num = i+block*j;

                for(int ii = start_x ; ii < end_x ; ii++)
                {
                    for(int jj = start_y ; jj < end_y ; jj++)
                    {
                        int index =data[jj][ii];
                        tmp2[num][index]++;
                    }
                }

                int average = width_block * height_block / 255;

                int LIMIT = 40 * average;
                int steal = 0;
                for(int k = 0 ; k < 256 ; k++)
                {
                    if(tmp2[num][k] > LIMIT){
                        steal += tmp2[num][k] - LIMIT;
                        tmp2[num][k] = LIMIT;
                    }
                }
                int bonus = steal/256;
                //hand out the steals averagely
                for(int k = 0 ; k < 256 ; k++)
                {
                    tmp2[num][k] += bonus;
                }

                for(int k = 0 ; k < 256 ; k++)
                {
                    if( k == 0)
                        C2[num][k] = 1.0f * tmp2[num][k] / total;
                    else
                        C2[num][k] = C2[num][k-1] + 1.0f * tmp2[num][k] / total;
                }
            }
        }
        //Calcular el valor del píxel transformado
        //Según la posición del píxel, elija un método de cálculo diferente
        for(int  i = 0 ; i < width; i++)
        {

            for(int j = 0 ; j < height; j++)
            {
                //four coners
                if(i <= width_block/2 && j <= height_block/2)
                {
                    int num = 0;
                    CLAHE_GO[j][i] = (int)(C2[num][CLAHE_GO[j][i]] * 255);
                }else if(i <= width_block/2 && j >= ((block-1)*height_block + height_block/2)){
                    int num = block*(block-1);
                    CLAHE_GO[j][i] = (int)(C2[num][CLAHE_GO[j][i]] * 255);
                }else if(i >= ((block-1)*width_block+width_block/2) && j <= height_block/2){
                    int num = block-1;
                    CLAHE_GO[j][i] = (int)(C2[num][CLAHE_GO[j][i]] * 255);
                }else if(i >= ((block-1)*width_block+width_block/2) && j >= ((block-1)*height_block + height_block/2)){
                    int num = block*block-1;
                    CLAHE_GO[j][i] = (int)(C2[num][CLAHE_GO[j][i]] * 255);
                }
                //four edges except coners
                else if( i <= width_block/2 )
                {

                    int num_i = 0;
                    int num_j = (j - height_block/2)/height_block;
                    int num1 = num_j*block + num_i;
                    int num2 = num1 + block;
                    float p =  (j - (num_j*height_block+height_block/2))/(1.0f*height_block);
                    float q = 1-p;
                    CLAHE_GO[j][i] = (int)((q*C2[num1][CLAHE_GO[j][i]]+ p*C2[num2][CLAHE_GO[j][i]])* 255);
                }else if( i >= ((block-1)*width_block+width_block/2)){

                    int num_i = block-1;
                    int num_j = (j - height_block/2)/height_block;
                    int num1 = num_j*block + num_i;
                    int num2 = num1 + block;
                    float p =  (j - (num_j*height_block+height_block/2))/(1.0f*height_block);
                    float q = 1-p;
                    CLAHE_GO[j][i] = (int)((q*C2[num1][CLAHE_GO[j][i]]+ p*C2[num2][CLAHE_GO[j][i]])* 255);
                }else if( j <= height_block/2 ){

                    int num_i = (i - width_block/2)/width_block;
                    int num_j = 0;
                    int num1 = num_j*block + num_i;
                    int num2 = num1 + 1;
                    float p =  (i - (num_i*width_block+width_block/2))/(1.0f*width_block);
                    float q = 1-p;
                    CLAHE_GO[j][i] = (int)((q*C2[num1][CLAHE_GO[j][i]]+ p*C2[num2][CLAHE_GO[j][i]])* 255);
                }else if( j >= ((block-1)*height_block + height_block/2) ){

                    int num_i = (i - width_block/2)/width_block;
                    int num_j = block-1;
                    int num1 = num_j*block + num_i;
                    int num2 = num1 + 1;
                    float p =  (i - (num_i*width_block+width_block/2))/(1.0f*width_block);
                    float q = 1-p;
                    CLAHE_GO[j][i] = (int)((q*C2[num1][CLAHE_GO[j][i]]+ p*C2[num2][CLAHE_GO[j][i]])* 255);
                }

                else{
                    int num_i = (i - width_block/2)/width_block;
                    int num_j = (j - height_block/2)/height_block;
                    int num1 = num_j*block + num_i;
                    int num2 = num1 + 1;
                    int num3 = num1 + block;
                    int num4 = num2 + block;
                    float u = (i - (num_i*width_block+width_block/2))/(1.0f*width_block);
                    float v = (j - (num_j*height_block+height_block/2))/(1.0f*height_block);
                    CLAHE_GO[j][i] = (int)((u*v*C2[num4][CLAHE_GO[j][i]] +
                        (1-v)*(1-u)*C2[num1][CLAHE_GO[j][i]] +
                        u*(1-v)*C2[num2][CLAHE_GO[j][i]] +
                        v*(1-u)*C2[num3][CLAHE_GO[j][i]]) * 255);
                }
                //smooth
                CLAHE_GO[j][i] = CLAHE_GO[j][i] + (CLAHE_GO[j][i] << 8) + (CLAHE_GO[j][i] << 16);
            }
        }
        dicomObj.saveData(CLAHE_GO,"/home/boris/Documentos/ClaheDicomH/csvC/MasaMicroC1_1.csv",",",false);
     return CLAHE_GO;
    }
