package semiglobalalignment;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.util.Scanner;
import java.lang.Math.*;
import UI.SemiGlobalAlignmentUI;
//import UI.hasil;

/**
 *
 * @author Kiky
 */

public class SemiGlobalAlignment {
    public static class fungsi{
    public static SemiGlobalAlignmentUI prim=new UI.SemiGlobalAlignmentUI();
    public static int[][] pam250;
    public static int[][] blosum62;
    public static int [][] dataTest;
    public static String dataQ; public static String[] DataQ;
    public static String dataD; public static String[] DataD;
    public static String[] sinkron={"ala","arg","asn","asp","cys","gln","glu","gly","his","ile","leu","lys","met","phe","pro","ser","thr","trp","tyr","val"};
    public static int gapInisialisasi,gapIterasi,baris, kolom, pilihangap, pilihanmatriks,max;
    public static enum Status{kiri, atas, diagonal};
    public static Status[][] status;
    public static String hasilAlignQ, hasilAlignD;
    public static boolean awalQ=false, akhirQ=false, awalD=false, akhirD=false; 
    public static void inisialisasi(){
        try{
           pam250=new int [20][20];
           blosum62=new int [20][20];
           FileInputStream in=new FileInputStream("pam250.txt");
           Scanner baca = new Scanner(in);
           BufferedReader line=new BufferedReader(new InputStreamReader(in));
           for (int i=0; i<20; i++){
               for (int j=0; j<20; j++)
                   pam250[i][j]=baca.nextInt();
               line.readLine();
           }
           FileInputStream in2=new FileInputStream("blosum62.txt");
           Scanner baca2 = new Scanner(in2);
           BufferedReader line2=new BufferedReader(new InputStreamReader(in2));
           for (int i=0; i<20; i++){
               for (int j=0; j<20; j++)
                   blosum62[i][j]=baca2.nextInt();
               line2.readLine();
           }
        }
        catch(Exception e){
            System.out.println(e);
        }
    }
    public static void matrikstest(){
        char[] dataQ2=dataQ.toCharArray();
        char[] dataD2=dataD.toCharArray();
        baris=(dataD2.length+1)/4+1;
        kolom=(dataQ2.length+1)/4+1;
        dataTest=new int [baris][kolom];
	if (awalQ==true){
	    for (int i=0; i<kolom; i++)
            dataTest[0][i]=0;
	}
	else{
	    for (int i=0; i<kolom; i++)
            dataTest[0][i]=i*(-1)*gapInisialisasi;
	}
	if (awalD==true){
	    for (int i=0; i<baris; i++)
            dataTest[i][0]=0;
	}
	else{
        for (int i=0; i<baris; i++)
            dataTest[i][0]=i*(-1)*gapInisialisasi;
	}
        DataQ=new String [kolom]; DataD=new String[baris];
        String temp="";
        DataQ[0]=null;DataD[0]=null;
        int i=1;
        System.out.print("\t");
        for (int k=0; k<dataQ2.length;k++){
                if (dataQ2[k]==' '){
                    DataQ[i]=temp;
                    System.out.print(DataQ[i]+"\t");
                    temp="";
                    i++;
                }
                else if(k==dataQ2.length-1){
                    temp += dataQ2[k];
                    DataQ[i]=temp;
                    System.out.print(DataQ[i]+"\n");
                    i=1;
                    temp="";
                }
                else
                    temp += dataQ2[k];
        }
        for (int k=0; k<dataD2.length;k++){
                if (dataD2[k]==' '){
                    DataD[i]=temp;
                    System.out.print(DataD[i]+"\n");
                    temp="";
                    i++;
                }
                else if(k==dataD2.length-1){
                    temp += dataD2[k];
                    DataD[i]=temp;
                    System.out.print(DataD[i]+"\n");
                    i=1;
                    temp="";
                }
                else
                    temp += dataD2[k];
        }
    }
    public static void algoritmaPAM250(){
        status=new Status[baris][kolom];
        status[0][0]=null;
        for (int i=1; i<baris; i++){
            status[i][0]=null;
            for (int j=1; j<kolom; j++){
                status[0][j]=null;
                int q=0,d=0,r;
                for (int k=0; k<20; k++){
                    if (DataQ[j].equals(sinkron[k]))
                        q=k;
                    if (DataD[i].equals(sinkron[k]))
                        d=k;
                }
                System.out.print(d+","+q+"\t");
                r=pam250[d][q];
                System.out.print(r+"\t");
                if (((dataTest[i-1][j-1]+r)>(dataTest[i][j-1]-gapIterasi))&&((dataTest[i-1][j-1]+r)>(dataTest[i-1][j]-gapIterasi))){
                    dataTest[i][j]=dataTest[i-1][j-1]+r;
                    status[i][j]=Status.diagonal;
                }
                else if(((dataTest[i][j-1]-gapIterasi)>(dataTest[i-1][j-1]+r))&&((dataTest[i][j-1]-gapIterasi)>(dataTest[i-1][j]-gapIterasi))){
                    dataTest[i][j]=dataTest[i][j-1]-gapIterasi;
                    status[i][j]=fungsi.Status.kiri;
                }
                else{
                    dataTest[i][j]=dataTest[i-1][j]-gapIterasi;
                    status[i][j]=fungsi.Status.atas;
                }
            }
            System.out.println();
        }
         for (int i=0; i<baris; i++){
            for (int j=0; j<kolom; j++){
                if(status[i][j]==fungsi.Status.diagonal)
                    System.out.print(dataTest[i][j]+""+status[i][j]+"\t");
                else
                    System.out.print(dataTest[i][j]+""+status[i][j]+"\t\t");
            }
            System.out.println();
         }   
    }
    public static void algoritmaBlosum62(){
        status=new Status[baris][kolom];
        for (int i=1; i<baris; i++){
            for (int j=1; j<kolom; j++){
                int q=0,d=0,r;
                for (int k=0; k<20; k++){
                    if (DataQ[j].equals(sinkron[k]))
                        q=k;
                    if (DataD[i].equals(sinkron[k]))
                        d=k;
                }
                System.out.print(d+","+q+"\t");
                r=blosum62[d][q];
                System.out.print(r+"\t");
                if (((dataTest[i-1][j-1]+r)>(dataTest[i][j-1]-gapIterasi))&&((dataTest[i-1][j-1]+r)>(dataTest[i-1][j]-gapIterasi))){
                    dataTest[i][j]=dataTest[i-1][j-1]+r;
                    status[i][j]=Status.diagonal;
                }
                else if(((dataTest[i][j-1]-gapIterasi)>(dataTest[i-1][j-1]+r))&&((dataTest[i][j-1]-gapIterasi)>(dataTest[i-1][j]-gapIterasi))){
                    dataTest[i][j]=dataTest[i][j-1]-gapIterasi;
                    status[i][j]=Status.kiri;
                }
                else{
                    dataTest[i][j]=dataTest[i-1][j]-gapIterasi;
                    status[i][j]=Status.atas;
                }
            }
            System.out.println();
        }
         for (int i=0; i<baris; i++){
            for (int j=0; j<kolom; j++){
                if(status[i][j]==Status.diagonal)
                    System.out.print(dataTest[i][j]+""+status[i][j]+"\t");
                else
                    System.out.print(dataTest[i][j]+""+status[i][j]+"\t\t");
            }
            System.out.println();
         }   
    }
    public static void alignment(){
	int d=0; int q=0;max=-100;
	boolean stop=false;
        hasilAlignD="";hasilAlignQ="";String gap="-----";String space="  ";
	if (akhirQ==true&&akhirD==true){
	    d=baris-1; q=kolom-1;
	}
	else if (akhirQ==true){
	    d=baris-1;
	    for (int j=1; j<kolom; j++){
		if (dataTest[d][j]>max){
		  q=j;
		}
	    }
	    for (int i=kolom-1; i>q; i--){
                hasilAlignQ=DataQ[i]+space+hasilAlignQ;
                hasilAlignD=gap+space+hasilAlignD;
            }
	}
	else if (akhirD==true){
	    q=kolom-1;
	    for (int j=1; j<baris; j++){
		if (dataTest[j][q]>max){
		    d=j;
		}
	    }
	    for (int i=baris-1; i>d; i--){
                hasilAlignD=DataD[i]+space+hasilAlignD;
                hasilAlignQ=gap+space+hasilAlignQ;
            }
	}
	else{
	    for (int j=1; j<kolom; j++){
		if (dataTest[baris-1][j]>max){
		    max=dataTest[baris-1][j];
		    d=baris-1;
		    q=j;
		}
	    }
	    for (int j=1; j<baris; j++){
		if (dataTest[j][kolom-1]>max){
		    max=dataTest[j][kolom-1];
		    d=j;
		    q=kolom-1;
		}
	    }
	    if (q!=kolom-1){
		for (int i=kolom-1; i>q; i--){
		    hasilAlignQ=DataQ[i].toLowerCase()+" "+hasilAlignQ;
		    hasilAlignD=space+" "+hasilAlignD;
		}
	    }
	    else if (d!=baris-1){
	        for (int i=baris-1; i>d; i--){
		    hasilAlignD=DataD[i].toLowerCase()+" "+hasilAlignD;
		    hasilAlignQ=space+" "+hasilAlignQ;
		}
	    }
	}
	max=dataTest[d][q];
        System.out.println(max+" "+status[d][q]+" "+d+","+q);
        /*if (q!=kolom-1){
            for (int i=kolom-1; i>q; i--){
                hasilAlignQ=DataQ[i]+space+hasilAlignQ;
                hasilAlignD=gap+space+hasilAlignD;
            }
        }
        else if (d!=baris-1){
            for (int i=baris-1; i>d; i--){
                hasilAlignD=DataD[i]+space+hasilAlignD;
                hasilAlignQ=gap+space+hasilAlignQ;
            }
        }*/
        while (stop==false){
            if (status[d][q]==Status.diagonal){
                hasilAlignQ=DataQ[q].toUpperCase()+space+hasilAlignQ;

                hasilAlignD=DataD[d].toUpperCase()+space+hasilAlignD;
                d=d-1;q=q-1;
            }
            else if (status[d][q]==Status.atas){
                hasilAlignQ=gap+space+hasilAlignQ;
                hasilAlignD=DataD[d].toUpperCase()+space+hasilAlignD;
                d=d-1;
            }
            else{
                hasilAlignQ=DataQ[q].toUpperCase()+space+hasilAlignQ;
                hasilAlignD=gap+space+hasilAlignD;
                q=q-1;
            }
            if(DataQ[q]==null){
                for (int i=d; i>0; i--){
                    hasilAlignQ=gap+space+hasilAlignQ;
                    hasilAlignD=DataD[i].toUpperCase()+space+hasilAlignD;
                }
                stop=true;
            }
            else if (DataD[d]==null){
                for (int i=q; i>0; i--){
                    hasilAlignQ=DataQ[i].toUpperCase()+space+hasilAlignQ;
                    hasilAlignD=gap+space+hasilAlignD;
                }
                stop=true;
            }
        }
	UI.hasil result=new UI.hasil();
	result.hasilD.setText(hasilAlignD);
	result.hasilQ.setText(hasilAlignQ);
	result.skor.setText(String.valueOf(max));
	awalD=false;awalQ=false;akhirQ=false;akhirD=false;
	result.setVisible(true);
        //System.out.println("q :\t"+hasilAlignQ+"\nd :\t"+hasilAlignD);
    }
}
    public static void main(String[] args) {
	fungsi.prim.setVisible(true);
    }
}
