package org.vmk.dep508.collections.library;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class LibraryImpl implements Library {

//    private Map<Book, BigInteger> library;

//    private Map<String, Map<String, BigInteger>> library;
    private Map<String, Map<Book, Integer>> library;

    public LibraryImpl(){
        this.library = new HashMap<String, Map<Book, Integer>>();
        this.library.put("library", new HashMap<Book, Integer>());
    }
    private Map<Book, Integer> getLibraryMap() {
        return this.library.get("library");
    }

    private void putBookToStr(Book book, String str){
        if (!library.containsKey(str)) {
            library.put(str, new HashMap<Book, Integer>());
        }

        if (library.get(str).containsKey(book)) {
            Integer tmp = library.get(str).get(book);
            tmp+= 1;
            library.get(str).put(book, tmp);
        } else {
            library.get(str).put(book, 1);
        }
    }

    private void getBookFromStr(Book book, String str){
        if (!library.containsKey(str)) {
            throw new BookNotFoundException(str + " is not exist");
        }
        if (library.get(str).containsKey(book) && (library.get(str).get(book) > 0) ) {
            Integer tmp = library.get(str).get(book);
            tmp = tmp - 1;
            library.get(str).put(book, tmp);
        } else {
            throw new BookNotFoundException("Book is not available");
//            library.get(str).put(book, 1);
        }
    }

    @Override
    public void addNewBook(Book book) {
       putBookToStr(book, "library");
    }

    @Override
    public void borrowBook(Book book, String student) {
        if (student.equals("library"))  {
            throw new BookNotFoundException(student +" - wrong student name");
        }
        getBookFromStr(book, "library");
        putBookToStr(book, student);
//        if (getLibraryMap().containsKey(book)) {
//            Integer tmp = getLibraryMap().get(book);
//            if (tmp.equals(0)) {
//                System.out.println("None books is left");
//                throw new BookNotFoundException("Book is not available");
//            } else {
//                tmp-=1;
//                System.out.println("Book was given to " + student);
//                System.out.println("Number of books now:" + tmp);
//                putBookToStr(book, student);
//            }
//        } else {
//            throw new BookNotFoundException("Book is not available");
//        }
    }

    @Override
    public void returnBook(Book book, String student) {
        if (student.equals("library"))  {
            throw new BookNotFoundException(student +" - wrong student name");
        }
        getBookFromStr(book, student);
        putBookToStr(book, "library");
//        if (getLibraryMap().containsKey(book)) {
//            Integer tmp = getLibraryMap().get(book);
//            if (tmp.equals(0)) {
//                System.out.println("None books is left");
//                throw new BookNotFoundException("Book is not available");
//            } else {
//                tmp-=1;
//                System.out.println("Book was given to " + student);
//                System.out.println("Number of books now:" + tmp);
//                putBookToStr(book, student);
//            }
//        } else {
//            throw new BookNotFoundException("Book is not available");
//        }
    }

    @Override
    public List<Book> findAvailableBooks() {
        List <Book> availableBooks = new ArrayList<Book>();

        for (Map.Entry<Book, Integer> entry : this.getLibraryMap().entrySet()){
            if (entry.getValue() >= 1) {
                availableBooks.add(entry.getKey());
            }
        }
        return availableBooks;

//        for (Map.Entry<String, Object> entry : map.entrySet()) {
//            String key = entry.getKey();
//            Object value = entry.getValue();
//            // ...
//        }
//        return null;
    }
}
