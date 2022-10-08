// Java code for thread creation by extending

// Task 1
class FirstTask extends Thread {
    public void run()
    {
        try {
            // Displaying the thread that is running
            int count = 1;
            while(true){
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

//Task 2
class SecondTask extends Thread {
    public void run()
    {
        try {
            // Displaying the thread that is running
            int count = 1;
            while(true){
                System.out.println("Thread " + Thread.currentThread().getId() + " is finishing Task 2 in loop " + count + " times");
                count ++;
                Thread.sleep(2000);
            }

        }
        catch (Exception e) {
            // Throwing an exception
            System.out.println("Exception is caught");
        }
    }
}

//Task 3
class ThirdTask extends Thread {
    public void run()
    {
        try {
            // Displaying the thread that is running
            int count = 1;
            while(true){
                System.out.println("Thread " + Thread.currentThread().getId() + " is finishing Task 3 in loop " + count + " times");
                count ++;
                Thread.sleep(3000);
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
    public static void main(String[] args)
    {
        int n = 3; // Number of threads
        for (int i = 0; i < n; i++) {
            if (i == 0){
                FirstTask firstTask = new FirstTask();
                firstTask.start();
            }
            else if (i == 1){
                SecondTask secondTask = new SecondTask();
                secondTask.start();
            }
            else if (i == 2){
                ThirdTask thirdTask = new ThirdTask();
                thirdTask.start();
            }
        }
    }
}