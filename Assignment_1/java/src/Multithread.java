// Java code for thread creation by extending

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

// Task 1
class FirstTask implements Runnable {
    Boolean threadDone;
    public FirstTask(){
    }

    @Override
    public void run()
    {
        try {
            // Displaying the thread that is running
            int count = 1;
            while(count <= 6){
                System.out.println("Thread " + Thread.currentThread().getId() + " is finishing Task 1 in loop " + count + " times");
                count ++;
                Thread.sleep(1000);
            }
        }
        catch (Exception e) {
            // Throwing an exception
            System.out.println("Exception is caught");
        }
    }
}


// Main Class
public class Multithread {
    public static void main(String[] args) throws InterruptedException {
        Boolean threadDone = false;
        int n = 3; // Number of threads
        FirstTask firstTask = new FirstTask();
        ExecutorService es = Executors.newCachedThreadPool();
        for (int i = 0; i < n; i++) {
            es.execute(new Runnable() {
                @Override
                public void run() {
                    try {
                        // Displaying the thread that is running
                        int count = 1;
                        while(count <= 6){
                            System.out.println("Thread " + Thread.currentThread().getId() + " is finishing Task 1 in loop " + count + " times");
                            count ++;
                            Thread.sleep(1000);
                        }
                    }
                    catch (Exception e) {
                        // Throwing an exception
                        System.out.println("Exception is caught");
                    }
                }
            });
        }
        es.shutdown();
        boolean finished = es.awaitTermination(1, TimeUnit.MINUTES);
        System.out.println("-----------------------------------------------------------------------------------------------------");
    }
}