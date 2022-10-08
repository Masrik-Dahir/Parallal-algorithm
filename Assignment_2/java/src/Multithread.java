// Java code for thread creation by extending

// Task 1
class FirstTask implements Runnable {
    int count;
    public FirstTask(int count) {
        this.count = count;
        // store parameter for later user
    }

    public void run() {
        System.out.println("fdsadfsadasd " + count);
    }
}

// Main Class
public class Multithread {
    public static void main(String[] args)
    {
        int n = 3; // Number of threads
        for (int i = 0; i < n; i++) {
            FirstTask firstTask = new FirstTask(i);
            new Thread(firstTask).start();
        }
    }
}