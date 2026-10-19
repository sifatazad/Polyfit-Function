#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define Abs(a) ((a)>=0?(a):-(a))

struct iteration {
    int *old_index;
    double *s;
    double *r;
};

static int compare_value(void const *a, void const *b) {
    struct iteration *it1= (struct iteration*) a;
    struct iteration *it2= (struct iteration*) b;

    if (*(it1->r) < *(it2->r))
        return 1;
    if (*(it1->r) > *(it2->r))
        return -1;
    return 0;
}

static int total_point(FILE *fp) {
    int n=0;
    char line[100];
    while (fgets(line,sizeof(line),fp)!=NULL) {
        n++;
    }
    rewind(fp);
    return n;
}

static double sum_x_pow(const double *x,int n,int p) {
    double s=0;
    for (int i=0; i<n; i++) {
        s+=pow(x[i],p*1.0);
    }
    return s;
}

static double sum_y_x_pow(const double *x,const double *y,const int n,const int p) {
    double s=0;
    for (int i=0; i<n; i++) {
        s+=pow(x[i],p*1.0)*y[i];
    }
    return s;
}

static void solve_equation_without_file(const int n,double **matrix,double *soln) {
    double *scaled_vec= malloc(sizeof(double)*n);
    if (!scaled_vec) {
        return;
    }
    for (int i=0; i<n; i++) {
        double mxm=0;
        for (int j=0; j<n; j++) {
            if (mxm<Abs(matrix[i][j]))
                mxm=Abs(matrix[i][j]);
        }
        scaled_vec[i]=mxm;
    }

    for (int i=0; i<n; i++) {
        int it=0;
        double *ratio=malloc(sizeof(double)*(n-i));
        if (!ratio) {
            free(scaled_vec);
            return;
        }
        double *new_scaled=malloc(sizeof(double)*(n-i));
        if (!new_scaled) {
            free(scaled_vec);
            free(ratio);
            return;
        }

        for (int j=i; j<n; j++) {
            new_scaled[it++]=scaled_vec[j];
        }
        it=0;
        for (int j=i; j<n; j++) {
            ratio[it]=Abs(matrix[j][i]/new_scaled[it]);
            it++;
        }

        int *old=malloc(sizeof(int)*(n-i));
        if (!old) {
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            return;
        }
        for (int k=i; k<n; k++) {
            old[k-i]=k;
        }

        struct iteration *iter=malloc(sizeof(struct iteration)*(n-i));
        if (!iter) {
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            free(old);
            return;
        }

        for (int j=0; j<n-i; j++) {
            iter[j].r=&ratio[j];
            iter[j].s=&new_scaled[j];
            iter[j].old_index=&old[j];
        }

        qsort(iter,n-i,sizeof(struct iteration),compare_value);

        int new=*iter[0].old_index;
        double *temp=malloc(sizeof(double)*(n+1));
        if (!temp) {
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            free(old);
            free(iter);
            return;
        }
        for (int x=0; x<n+1; x++) {
            temp[x]=matrix[i][x];
        }
        for (int x=0; x<n+1; x++) {
            matrix[i][x]=matrix[new][x];
        }
        for (int x=0; x<n+1; x++) {
            matrix[new][x]=temp[x];
        }
        double val=matrix[i][i];
        if (Abs(val)<1e-9) {
            printf("Singular Matrix");
            free(scaled_vec);
            free(ratio);
            free(new_scaled);
            free(old);
            free(iter);
            free(temp);
            return;
        }
        for (int r=0; r<n+1; r++) {
            matrix[i][r]/=val;
        }

        double tmp=scaled_vec[new];
        scaled_vec[new]=scaled_vec[i];
        scaled_vec[i]=tmp;

        for (int k=0; k<n; k++) {
            if (i==k)
                continue;
            double factor= -matrix[k][i]/matrix[i][i];
            for (int col=0; col<n+1; col++) {
                matrix[k][col]+=factor*matrix[i][col];
            }
        }

        free(ratio);
        free(iter);
        free(old);
        free(temp);
        free(new_scaled);
    }
    for (int i=0; i<n; i++) {
        soln[i]=matrix[i][n];
    }
    // free_matrix(matrix,n);
    free(scaled_vec);
}

static void polyfit(FILE *fp,double *soln,int n,int fit_rank) {
    double **matrix= malloc(sizeof(double*)*(fit_rank+1));
    if (!matrix) {
        return;
    }
    for (int i=0; i<=fit_rank; i++) {
        matrix[i]=malloc(sizeof(double)*(fit_rank+2));
        if (!matrix[i]) {
            free(matrix);
            return;
        }
    }

    double *x,*y;
    x=malloc(sizeof(double)*n);
    if (!x) {
        for (int i=0; i<=fit_rank; i++) {
            free(matrix[i]);
        }
        free(matrix);
        return;
    }
    y=malloc(sizeof(double)*n);
    if (!y) {
        for (int i=0; i<=fit_rank; i++) {
            free(matrix[i]);
        }
        free(matrix);
        free(x);
        return;
    }

    for (int i=0; i<n; i++) {
        fscanf(fp,"%lf\t%lf",&x[i],&y[i]);
    }

    for (int i=0; i<=fit_rank; i++) {
        for (int j=0; j<=fit_rank; j++) {
            matrix[i][j]=sum_x_pow(x,n,i+j);
        }
    }

    for (int i=0; i<=fit_rank; i++) {
        matrix[i][fit_rank+1]=sum_y_x_pow(x,y,n,i);
    }

    FILE *fo=fopen("C:/Users/Rifat/OneDrive/Documents/polyfit/output.txt","w");
    if (!fo) {
        for (int i=0; i<=fit_rank; i++) {
            free(matrix[i]);
        }
        free(matrix);
        free(x);
        free(y);
        return;
    }

    for (int i=0; i<=fit_rank; i++) {
        for (int j=0; j<=fit_rank+1; j++) {
            fprintf(fo,"%.3lf\t",matrix[i][j]);
        }
        fprintf(fo,"\n");
    }

    solve_equation_without_file(fit_rank+1,matrix,soln);

    // freeing here
    for (int i=0; i<=fit_rank; i++) {
        free(matrix[i]);
    }
    free(matrix);
    free(x);
    free(y);
    fclose(fo);
}

static void polyfit_without_file(const double *x, const double *y, double *soln,
    int n, int fit_rank) {
    double **matrix= malloc(sizeof(double*)*(fit_rank+1));
    if (!matrix) {
        return;
    }
    for (int i=0; i<=fit_rank; i++) {
        matrix[i]=malloc(sizeof(double)*(fit_rank+2));
        if (!matrix[i]) {
            free(matrix);
            return;
        }
    }

    for (int i=0; i<=fit_rank; i++) {
        for (int j=0; j<=fit_rank; j++) {
            matrix[i][j]=sum_x_pow(x,n,i+j);
        }
    }

    for (int i=0; i<=fit_rank; i++) {
        matrix[i][fit_rank+1]=sum_y_x_pow(x,y,n,i);
    }

    // FILE *fo=fopen("C:/Users/Rifat/OneDrive/Documents/polyfit/output.txt","w");
    // if (!fo) {
    //     for (int i=0; i<=fit_rank; i++) {
    //         free(matrix[i]);
    //     }
    //     free(matrix);
    //     return;
    // }
    //
    // for (int i=0; i<=fit_rank; i++) {
    //     for (int j=0; j<=fit_rank+1; j++) {
    //         fprintf(fo,"%.3lf\t",matrix[i][j]);
    //     }
    //     fprintf(fo,"\n");
    // }

    solve_equation_without_file(fit_rank+1,matrix,soln);

    // freeing here
    for (int i=0; i<=fit_rank; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

int main() {
    FILE *fp=fopen("C:/Users/Rifat/OneDrive/Documents/polyfit/input.txt","r");
    if (!fp) {
        return 1;
    }
    int tot= total_point(fp);

    int n=5;
    double *co_eff= malloc(sizeof(double)*(n+1));
    if (!co_eff) {
        fclose(fp);
        return 1;
    }
    polyfit(fp,co_eff,tot,n);

    for (int i=0; i<=n; i++) {
        printf("%.6lf\n",co_eff[i]);
    }

    fclose(fp);
    free(co_eff);
}