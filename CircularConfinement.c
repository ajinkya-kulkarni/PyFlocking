#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h> 

#define tf 14000      
#define dt 0.01         
#define datafreq 100    
#define beta 0.01       
#define gamma 0.0001    

#define totalpart 7289 
#define core 6120            
#define ntotalin 6120
#define xoutw 1169      

#define cellfreq 1     
int sstate=0,sstate1=0,sstate2=0,sstate3=0,sstate4=0,sstate5=0,steadytimesteps=0; 
int startsstate=0,startsstate1=0,startsstate2=0,startsstate3=0,startsstate4=0,startsstate5=0,packet=0;
#define ga 0         
#define rho 100     
#define hsquare 61.0351    
#define h 7.8125           
#define cell_length 2.5  
float dm = 0.500;    
float mass = 59.999047; 
float cv;            
float video_time;   
float limit_overlap=0;  
float omega=0;     
float pi=3.14159265358979323846264338327; 

float r_cavity1=22.57; 
float x_c=40; 	   
float y_c=40;      

#define length 100  
float domain_length;
int cellrow=length/cell_length; 
int cellcol=length/cell_length; 
int nbrpart=200;     
int celltotal;    
int cell_totalpart[1600][200]; 


#define sigma 305.58 
float V_surr = hsquare;  
int timestep;
void force_calc(float fx[],float fy[],float x[],float y[],float vx[],float vy[], float k_stiff[], float cv, float Fc, float Fd, float Fsp, float Fnet);
int cj=0,ck=0,cf=0;
int xloc,yloc,cellid;   
int n_influ=1471;
float wij[1471];           
float ds;                 

  float vx[totalpart],vy[totalpart];  
  float vonep[totalpart];             
  float ronep[totalpart];             
  float vavg_onet,vsum_onet=0.0;      
  float vronet,vthonet,vrsum_onet=0.0,vthsum_onet=0.0,vrallt=0.0,vthallt=0.0,vcore_onet=0.0,vfree_onet=0.0,vcore_allt=0.0;
  float vronet1,vthonet1,vrsum_onet1=0.0,vthsum_onet1=0.0,vcore_onet1=0.0,vfree_onet1=0.0;
  float vronet2,vthonet2,vrsum_onet2=0.0,vthsum_onet2=0.0,vcore_onet2=0.0,vfree_onet2=0.0;
  float vronet3,vthonet3,vrsum_onet3=0.0,vthsum_onet3=0.0,vcore_onet3=0.0,vfree_onet3=0.0;
  float vronet4,vthonet4,vrsum_onet4=0.0,vthsum_onet4=0.0,vcore_onet4=0.0,vfree_onet4=0.0;
  float vrallt1=0.0,vrallt2=0.0,vrallt3=0.0,vrallt4=0.0;
  float vthallt1=0.0,vthallt2=0.0,vthallt3=0.0,vthallt4=0.0;
  float vcore_allt1=0.0,vcore_allt2=0.0,vcore_allt3=0.0,vcore_allt4=0.0;
  float vavg1,vavg2,vavg3,vavg4;
  float vavg,vsum_allt=0.0;           
  float THETA,vronep[totalpart],vthonep[totalpart];
  float x[totalpart],y[totalpart];   
  float xi[totalpart],yi[totalpart];  
  float vxi[totalpart],vyi[totalpart];
  float k_stiff[totalpart];           
  float time1,dx,dy;
  float fxi[totalpart],fyi[totalpart];
  float fx[totalpart],fy[totalpart]; 
  int i=0,timestep,k,a,tstepfinal;
  float d=0;
  float di=0.005;

float Rring=11.285; 
float halfdeltar=0.56425; 
float Rringm=10.72075;
float Rringp=11.84925;
float ringcount=0.0,ringcountonet=0.0,ringcountallt=0.0,ringcountCUM=0.0,ringcountAVG=0.0,ringcountalltsq=0.0,ringcountalltsqCUM=0.0,sigmasq=0.0,nsqEnsemble=0.0,sqEnsemblen=0.0;

  int i,j,k,id;
  float dx,dy,d,f;
  float k_eq; 
  float vel_mag, v_i = 0.0, v_j = 0.0, vm_i = 0.0, vm_j = 0.0;
  float fcx[totalpart],fcy[totalpart],Fc;
  float fdx[totalpart],fdy[totalpart],Fd,fspx[totalpart],fspy[totalpart],Fsp,Fnet;
  int cell[9],cal[500];
  int count;
  float vxcp1,vcp2,vycp1,prod;
  float m_part, m_fluid;
  int xcoord,ycoord,nbrvect;
  float xp,yp,theta;
  float overlap;

#define startsstate  4000 
      
#define startsstate1 4000
#define startsstate2 6500
#define startsstate3 9000
#define startsstate4 11500

#define	packet  2500
#define sstate1 2500
#define sstate2 2500
#define sstate3 2500
#define sstate4 2500


/********************* start of main() function   ************************/

int main()
{
  FILE *LOC, *VEL;             
  FILE *VELsstate, *VELavg;
  FILE *Vonet;  
  FILE *fcompalltime; 	      
  FILE *fnetalltime; 	      
  FILE *iLOC,*iVEL;             
  FILE *cvin;                
  FILE *VELcore;             
  FILE *VELavgcore;
FILE *Nringonet, *Nringallt, *NringAVG, *NsqringAVG, *GNF;
 
  for(i=0;i<n_influ;i++)
    {
      wij[i]=exp(-(d*d)/h); 
      d=d+di;
    }

  iLOC=fopen("init_positions.txt","r"); 
  iVEL=fopen("init_velocities.txt","r"); 
  cvin=fopen("cv.txt","r");              

  for(i=0;i<totalpart;i++)
    {
      if(fscanf(iLOC,"%d %f %f %f\n",&i,&x[i],&y[i],&k_stiff[i])){};
      if(fscanf(iVEL,"%d %f %f\n",&i,&vx[i],&vy[i])){};
    }

  if(fscanf(cvin,"%f",&cv)){};  

  fclose(iLOC);
  fclose(iVEL);
  fclose(cvin);

	dx=length/nbrpart;  
	dy=length/nbrpart;  
	celltotal=cellrow*cellcol;

	force_calc(fx,fy,x,y,vx,vy,k_stiff,cv,Fc,Fd,Fsp,Fnet);

	LOC=fopen("data.txt","w");         
	VEL=fopen("velocity.txt","w");     
	VELsstate=fopen("Vpolar_sstate_GAMMA.txt","w"); 
        Vonet=fopen("vel_rtheta_allt.txt","w");
        VELavg=fopen("v_AVG_rtheta.txt","w"); 
        VELcore=fopen("vel_xy_all_sstate.txt","w");
        VELavgcore=fopen("v_AVG_xy.txt","w"); 
        fcompalltime=fopen("fcomp_alltime.txt","w");
	fnetalltime=fopen("fnet_alltime.txt","w");
Nringonet=fopen("ring_count_onetime.txt","w");
Nringallt=fopen("ring_count_alltime.txt","w");
NringAVG=fopen("ring_count_AVG.txt","w");
NsqringAVG=fopen("square_of_ring_count_AVG.txt","w");
GNF=fopen("GNF.txt","w");

if(omega==0)
	{
	for(i=ntotalin;i<ntotalin+xoutw;i++) 
		{
			vx[i] = 0;
			vy[i] = 0;
		}
	}

	tstepfinal = tf/dt + 1;  
	sstate=tf-startsstate;
        steadytimesteps = sstate/dt;  
for(timestep=0;timestep<tstepfinal;timestep++)
      { 
      time1=(timestep+1)*dt;
      cf=timestep/cellfreq;
      cj=timestep;
      a=timestep/datafreq;

      if(timestep==a*datafreq)
	{
	  video_time = timestep*dt;

          if (omega==0)
          {
	   printf("Core particles are self propelled along predetermined motive direction(V0) of strength GAMMA: \n");	   
	   printf("   radius of cavity = %0.2fm\n",r_cavity1);
           printf("Particles simulated = %d\n",core);
	   domain_length=length;
           printf("cell domain length  = %0.2fm\n",domain_length);
           printf("              GAMMA = %0.4f\n",gamma);
	   printf("               BETA = %0.4f\n",beta);        
           printf("Forward trial at CV = %0.2f\n",cv);
           printf("Total Flow time     = %d secs\n",tf);    
           printf("Video generated     = %0.6f secs\n",video_time);
           printf("\n");
	  }

	  vrsum_onet=0.0;
          vthsum_onet=0.0;
          vcore_onet=0.0;
          ringcountonet=0.0;
          ringcountallt=0.0;
	  for(k=0;k<totalpart;k++)
	    {
	      vonep[k]=sqrt(((vx[k]*vx[k])+(vy[k]*vy[k]))); 
 	      ronep[k]=sqrt((((x[k]-x_c)*(x[k]-x_c))+((y[k]-y_c)*(y[k]-y_c)))); 

 	      fprintf(LOC,"%f\t%f\t%f\t%f\t%f\t%d\n",time1-dt,x[k],y[k],dm/2,ronep[k],k+1);

	      fprintf(VEL,"%f\t%f\t%f\t%f\t%d\n",time1-dt,vx[k],vy[k],vonep[k],k+1);

                 if(timestep>=startsstate/dt) 
                 {
                 if(k<ntotalin)
                  {
					yp=y[k]-y_c;
					xp=x[k]-x_c;
					THETA = atan2(yp, xp);
					if (THETA< 0)
                                           {  
                                            THETA = THETA + 2*pi;
                                           }					
 				        vronep[k] = vy[k]*sin(THETA) + vx[k]*cos(THETA);
                                        vthonep[k] = vy[k]*cos(THETA) - vx[k]*sin(THETA);

 	           fprintf(VELsstate,"%d\t%d\t%f\t%f\t%f\n",cj,k+1,ronep[k],vronep[k],vthonep[k]);
		   vrsum_onet = vrsum_onet + vronep[k];
		   vthsum_onet = vthsum_onet + vthonep[k];    		   
		   if(ronep[k]>Rringm && ronep[k]<Rringp)
		   ringcountonet=ringcountonet+((ronep[k])/ronep[k]);
                   fprintf(Nringonet,"%f\n",ringcountonet);
                   
	          }
		}
	      }

              if(timestep>=startsstate/dt) 
                {
                vronet = (vrsum_onet*dt*datafreq)/(core);
                vthonet= (vthsum_onet*dt*datafreq)/(core);
		               
		vrallt=vrallt+vronet;
		vthallt=vthallt+vthonet;
                fprintf(Vonet,"%d\t%f\t%f\t%f\t%f\n",cj,vronet,vthonet,vrallt,vthallt);	                
                ringcountallt=(ringcountonet+ringcountallt);
                ringcountalltsq=(ringcountallt*ringcountallt);
                ringcountCUM=ringcountCUM + ringcountallt;
                ringcountalltsqCUM=(ringcountalltsqCUM + ringcountalltsq);
                fprintf(Nringallt,"%f\n",ringcountallt);
                }


             if(timestep>=startsstate/dt) 
                {
                 for(k=0;k<ntotalin;k++)
	            {
       	             vonep[k]=sqrt(((vx[k]*vx[k])+(vy[k]*vy[k]))); 
 	             fprintf(VELcore,"%f\t%d\t%f\t%f\t%f\n",time1-dt,k+1,vx[k],vy[k],vonep[k]);
 	             vcore_onet = vcore_onet + vonep[k];
	            } 
	            vfree_onet = (vcore_onet*dt*datafreq)/(core);
	            vcore_allt = vcore_allt+vfree_onet;
                }

	  fprintf(LOC,"\n");
	  fprintf(VEL,"\n");
	  fprintf(VELsstate,"\n");
	  fprintf(Vonet,"\n");
          fprintf(VELcore,"\n");
          fprintf(Nringonet,"\n");
          fprintf(Nringallt,"\n");          
	}  
      for(i=0;i<ntotalin;i++) 
         {

     	  x[i]=x[i]+(vx[i]*dt)+(((dt*dt)*(fx[i]/mass))/2); 
      	  y[i]=y[i]+(vy[i]*dt)+(((dt*dt)*(fy[i]/mass))/2);

	  vxi[i]=vx[i]; 
	  vyi[i]=vy[i]; 

	  vx[i]=vx[i]+((fx[i]/mass)*dt); 
	  vy[i]=vy[i]+((fy[i]/mass)*dt); 

	  fxi[i]=fx[i];
	  fyi[i]=fy[i];

         }
	force_calc(fx,fy,x,y,vx,vy,k_stiff,cv,Fc,Fd,Fsp,Fnet);

            {
              for(i=2000;i<2001;i++) 
                 {
                  Fc=sqrt(((fcx[i]*fcx[i]) + (fcy[i]*fcy[i])));
		  
                  Fd=sqrt(((fdx[i]*fdx[i]) + (fdy[i]*fdy[i])));

                  Fsp=sqrt(((fspx[i]*fspx[i]) + (fspy[i]*fspy[i])));

		  Fnet=sqrt(((fx[i]*fx[i]) + (fy[i]*fy[i])));

		  fprintf(fcompalltime,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%d",cj,fcx[i],fcy[i],fdx[i],fdy[i],fspx[i],fspy[i],i);

		  fprintf(fnetalltime,"%d\t%f\t%f\t%f\t%f\t%d",cj,Fc,Fd,Fsp,Fnet,i);
	         }

             }

                 fprintf(fcompalltime,"\n");
                 fprintf(fnetalltime,"\n");

      for(i=0;i<ntotalin;i++)
	 {
	  vx[i]=vxi[i]+((fxi[i]/mass)+(fx[i]/mass))*dt*0.5; 
	  
	  vy[i]=vyi[i]+((fyi[i]/mass)+(fy[i]/mass))*dt*0.5; 
	  
	 }

     }
     
                 { 
                  vrallt=(vrallt)/sstate;
                  vthallt=(vthallt)/sstate;
                  vavg=sqrt(((vrallt*vrallt) + (vthallt*vthallt)));              
                  fprintf(VELavg,"%f\n",vavg);  		
  		 } 

               fprintf(VELavg,"\n");


                 { 
                  vcore_allt=(vcore_allt)/(sstate);                                          
                  fprintf(VELavgcore,"%f\n",vcore_allt);  		
  		 } 
               fprintf(VELavgcore,"\n");                  

{ 
ringcountAVG = ((ringcountCUM)/(sstate));
nsqEnsemble = ((ringcountalltsqCUM)/(sstate)); 
sqEnsemblen = (ringcountAVG*ringcountAVG);
sigmasq     = ((nsqEnsemble - sqEnsemblen)/(ringcountAVG));
}
fprintf(NringAVG,"%f\n",ringcountAVG);
fprintf(NsqringAVG,"%f\n",nsqEnsemble);
fprintf(GNF,"%f\n",sigmasq);

  	fclose(LOC);
  	fclose(VEL);
	fclose(VELsstate);
        fclose(Vonet);
        fclose(VELavg);
        fclose(fcompalltime);
        fclose(fnetalltime);
        fclose(VELcore);
        fclose(VELavgcore);
        fclose(Nringonet);
        fclose(Nringallt);
        fclose(NringAVG);
        fclose(NsqringAVG);
        fclose(GNF);
}
/**************************  END of main() function   ********************/


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    FORCE CALCULATION function   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

void force_calc(float fx[],float fy[],float x[],float y[],float vx[],float vy[],float k_stiff[],float cv,float Fc,float Fd,float Fsp,float Fnet)

if(cf==cj*cellfreq)
    {
      for(i=0;i<celltotal;i++) 

	{
	  cell_totalpart[i][1]=0; 

	  for(k=2;k<nbrpart;k++) 
	    {
	      cell_totalpart[i][k]=-1; 
	    }
	}

      for(i=0;i<totalpart;i++) 
	{
	    xloc=x[i]/cell_length;
	    yloc=y[i]/cell_length;

	  cellid=((yloc+1)*cellcol)+(xloc+1); 
	  
	  cell_totalpart[cellid][1]=cell_totalpart[cellid][1]+1; 

	  xloc=cell_totalpart[cellid][1];
	  cell_totalpart[cellid][xloc+1]=i; 
	}
    }

  for(i=0;i<totalpart;i++) 
    {
      fcx[i]=0;
      fcy[i]=-mass*ga;
         
        xcoord=x[i]/cell_length;
        ycoord=y[i]/cell_length;
      
      id=((ycoord+1)*cellcol)+(xcoord+1); 

      cell[0]=id+cellrow+1;
      cell[1]=id+cellrow;
      cell[2]=id+cellrow-1;
      cell[3]=id+1;
      cell[4]=id;
      cell[5]=id-1;
      cell[6]=id-(cellrow-1);
      cell[7]=id-cellrow;
      cell[8]=id-(cellrow+1);

      count=0; 

      for(k=0;k<9;k++) 
	{
          if(cell_totalpart[cell[k]][1]>0) 
	    {
	      for(j=1;j<=cell_totalpart[cell[k]][1];j++)
		{
		  cal[j-1+count]=cell_totalpart[cell[k]][j+1];
		}
	      count=count+cell_totalpart[cell[k]][1]; 
	    }
	}

/*##################################  1 PARTICLE PARTICLE FORCE   ######################################*/

      vxcp1=0;
      vcp2=0;
      vycp1=0;
      m_part=0;

      for(j=0;j<count;j++) 
	{
   	   if(i!=cal[j]) 
	    {
	      dx=x[cal[j]]-x[i]; 

	      dy=y[cal[j]]-y[i];

   	       d=(sqrt(((dx*dx)+(dy*dy))));

	     ds = dm;

	overlap = (ds-d);	

 if(d<0.2*ds) 
	{
	  d=0.1*ds; 
	}

 if(overlap>limit_overlap)
	{	
       	  k_eq = (k_stiff[i]*k_stiff[cal[j]])/(k_stiff[i] + k_stiff[cal[j]]);

	  f=(-k_eq*(ds-d));

	 fcx[i]=fcx[i]+(f*(dx/d));  

	 fcy[i]=fcy[i]+(f*(dy/d));  

	 Fc=sqrt(((fcx[i]*fcx[i]) + (fcy[i]*fcy[i])));
	}

nbrvect=1+n_influ*(d/h);
prod=wij[nbrvect];
vcp2=vcp2 + prod*mass; 
m_part = m_part + mass;
vxcp1=vxcp1 + prod*(mass*vx[cal[j]]); 
vycp1=vycp1 + prod*(mass*vy[cal[j]]); 

	  }
      }


/*####################################   2 DRAG FORCE   ########################################*/

m_fluid = rho*(V_surr - m_part/sigma);

      if(vcp2==0)
	{
	  fdx[i]=cv*dm*(0 - vx[i]); 
	  fdy[i]=cv*dm*(0 - vy[i]); 
	}

      if(vcp2!=0)
	{
	  fdx[i]=cv*dm*((m_part/(m_part + m_fluid))*(vxcp1/vcp2) - vx[i]); 
	  fdy[i]=cv*dm*((m_part/(m_part + m_fluid))*(vycp1/vcp2) - vy[i]); 
	}

          Fd=sqrt(((fdx[i]*fdx[i]) + (fdy[i]*fdy[i])));

/*##############################   3 SELF PROPELLED FORCE   ################################*/

				  vel_mag = sqrt(((vx[i]*vx[i]) + (vy[i]*vy[i])));

       				  if (vel_mag > 1e-30)
          			  {
            				v_i = vx[i]/sqrt(((vx[i]*vx[i]) + (vy[i]*vy[i])));
            				v_j = vy[i]/sqrt(((vx[i]*vx[i]) + (vy[i]*vy[i])));
            			  }

				
          			  {
					yp=y[i]-y_c;
					xp=x[i]-x_c;

					theta = atan2(yp, xp);
					if (theta < 0) theta = theta + 2*pi;
					vm_i = cos(theta+pi/2);
					vm_j = sin(theta+pi/2);
					
    				  }
					fspx[i]=mass*(beta*v_i + gamma*vm_i); 
					fspy[i]=mass*(beta*v_j + gamma*vm_j); 
                                        Fsp=sqrt(((fspx[i]*fspx[i]) + (fspy[i]*fspy[i])));

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  final  FORCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

	         fx[i]=fcx[i]+ fdx[i] + fspx[i]; 
		 fy[i]=fcy[i]+ fdy[i] + fspy[i]; 
                  Fnet=sqrt(((fx[i]*fx[i]) + (fy[i]*fy[i])));
     }
}

/*%%%%%%%%%%%%%%%%% END of  FORCE CALCULATION function  %%%%%%%%%%%%%%*/

