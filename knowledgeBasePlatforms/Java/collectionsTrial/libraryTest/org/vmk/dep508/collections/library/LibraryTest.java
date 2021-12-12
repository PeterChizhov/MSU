package org.vmk.dep508.collections.library;

import org.junit.After;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import java.util.List;

import static org.junit.Assert.*;
import static org.hamcrest.CoreMatchers.*;

public class LibraryTest {

    Library libr;
//    String product4thEdition = "Thinking Java 4th edition";
//    String product3rdEdition = "Thinking Java 3rd edition";
//    BigDecimal defaultPrice = BigDecimal.TEN;
    Book book1 = new Book(1, "Hello");
    Book book2 = new Book(2, "World");
    Book book3 = new Book(3, "!");
    Book book4 = new Book(3, "!!!!!");

    /*Вызыватется перед вызовом каждого метода помеченного аннотацией @Test*/
    @Before
    public void setup(){
        libr = new LibraryImpl();
        libr.addNewBook(book1);
        libr.addNewBook(book1);
        libr.addNewBook(book1);
        libr.addNewBook(book2);
        libr.addNewBook(book3);
    }

    /*Вызыватется после вызова каждого метода помеченного аннотацией @Test*/
    @After
    public void clear() {
        libr = null;
    }

    @Test
    public void checkCollection() {
        List<Book> bookLst = libr.findAvailableBooks();

        assertThat(bookLst, hasItems(book1, book2, book3));

        Book tmp = book4;
        assertThat(bookLst, not(hasItem(tmp)));


        libr.borrowBook(book1, "Kirill");
        libr.borrowBook(book1, "Ramil");
        libr.borrowBook(book1, "Ramil");
        List<Book> bookLst1 = libr.findAvailableBooks();

        assertThat(bookLst1, hasItems(book2, book3));
        assertThat(bookLst1, not(hasItems(book1)));

    }

    @Test(expected = BookNotFoundException.class)
    public void addNewBook() throws Exception{
        List<Book> bookLst = libr.findAvailableBooks();
        assertThat(bookLst, hasItems(book1, book2, book3));
        Book tmp = book4;
        assertThat(bookLst, not(hasItem(tmp)));
        libr.borrowBook(book3, "library");
    }

    @Test(expected = BookNotFoundException.class)
    public void borrowBook() {
        libr.borrowBook(book4, "Ramil");
    }

    @Test(expected = BookNotFoundException.class)
    public void returnBook() {
        libr.borrowBook(book3, "Ramil");
        libr.borrowBook(book2, "Ramil");
        libr.borrowBook(book1, "Ramil");

        List<Book> bookLst = libr.findAvailableBooks();
        assertThat(bookLst, hasItems(book1));
        assertThat(bookLst, not(hasItems(book2, book3)));

        libr.returnBook(book1, "Ramil");
        libr.returnBook(book2, "Ramil");
        libr.returnBook(book3, "Ramil");

        List<Book> bookLst1 = libr.findAvailableBooks();
        assertThat(bookLst1, hasItems(book1, book2, book3));

        libr.returnBook(book1, "Ramil");
    }

}