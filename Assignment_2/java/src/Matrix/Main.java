package Matrix;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

interface MatMath{
    void multiply(int[][] a, int[][] b, int[][] c) throws InterruptedException;     //multiply A and B into C
    void add(int[][] a, int[][] b, int[][] c) throws InterruptedException;          //add A and B into C
    void print(int[][] a) throws InterruptedException;                              //pretty print A
}
class RowColAddExecutable implements Runnable{

    private CountDownLatch doneSignal;
    int[][] first, second, result;
    int row, col;

    public RowColAddExecutable(CountDownLatch doneSignal, int[][] first, int[][] second, int[][] result, int row, int col){
        this.doneSignal = doneSignal;
        this.first = first;
        this.second = second;
        this.result = result;
        this.row = row;
        this.col = col;
    }

    public void run(){
        doWork();
        doneSignal.countDown();
    }

    private void doWork(){
        result[row][col] += first[row][col] + second[row][col];
    }
}
class RowColProdExecutable implements Runnable{

    private CountDownLatch doneSignal;
    int[][] first, second, result;
    int row, col, size;

    public RowColProdExecutable(CountDownLatch doneSignal, int[][] first, int[][] second, int[][] result, int row, int col, int size){
        this.doneSignal = doneSignal;
        this.first = first;
        this.second = second;
        this.result = result;
        this.row = row;
        this.col = col;
        this.size = size;
    }

    public void run(){
        doWork();
        doneSignal.countDown();
    }

    private void doWork(){
        for(int k = 0; k < size; k++){
            result[row][col] += first[row][k] * second[k][col];
        }
    }
}
class ParallelMatMath implements MatMath {

    public ParallelMatMath(){}

    public void multiply(int[][] a, int[][] b, int[][] c) throws InterruptedException{
        CountDownLatch doneSignal = new CountDownLatch(a.length * b[0].length);
        /*Executor e = new Executor() {
            @Override
            public void execute(Runnable command) { new Thread(command).start(); }
        };*/
        ExecutorService e = Executors.newFixedThreadPool(a.length * b[0].length);
        for(int i = 0; i < a.length; i++){
            for(int j = 0; j < b[0].length; j++) {
                e.execute(new RowColProdExecutable(doneSignal, a, b, c, i, j, a.length));
            }
        }
        doneSignal.await();
    }

    public void add(int[][] a, int[][] b, int[][] c)throws InterruptedException{
        CountDownLatch doneSignal = new CountDownLatch(a.length * a[0].length);
        /*Executor e = new Executor() {
            @Override
            public void execute(Runnable command) { new Thread(command).start(); }
        };*/
        ExecutorService e = Executors.newFixedThreadPool(a.length * a[0].length);
        for(int i = 0; i < a.length; i++){
            for(int j = 0; j < b[i].length; j++) {
                e.execute(new RowColAddExecutable(doneSignal, a, b, c, i, j));
            }
        }
        doneSignal.await();
    }

    public void print(int[][] a)throws InterruptedException{
        for(int i = 0; i < a.length; i++){
            for(int j = 0; j < a[i].length; j++){
                System.out.printf("%5d", a[i][j]);
            }
            System.out.println();
        }
    }
}

public class Main {

    public static int[][] fillMatrix(int[][] m){
        for(int i = 0; i < m.length; i++){
            for(int j = 0; j < m[i].length; j++){
                m[i][j] = (int)(Math.random()*9);
            }
        }
        return m;
    }

    public static void main(String[] args) {

        int[][] A, B, C, D, r, s, t;

        // comment/uncomment for different implemenations
        MatMath p = new ParallelMatMath();
        //Matrix.MatMath p = new ParaStreamMatMath();

        //code for ABCD
        A = fillMatrix(new int[4][4]);
        B = fillMatrix(new int[4][4]);
        t = new int[A.length][B[0].length];


        try{

            System.out.println("Matrix A(" + A.length + " x " + A[0].length + "):");
            p.print(A);
            System.out.println();

            System.out.println("Matrix B(" + B.length + " x " + B[0].length + "):");
            p.print(B);
            System.out.println();


            //multiple test
            System.out.println("Multiply A * B = t");
            p.multiply(A, B, t);

            System.out.println("Matrix t(" + t.length + " x " + t[0].length + "):");
            p.print(t);
            System.out.println();


        }catch (InterruptedException ex){
            System.out.print(ex.getStackTrace().toString());
        }

    }
}