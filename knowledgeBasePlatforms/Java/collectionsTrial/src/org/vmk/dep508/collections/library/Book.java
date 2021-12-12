package org.vmk.dep508.collections.library;

import java.util.Objects;

public class Book {
    final int id;
    final String title;

    public Book(int id, String title) {
        this.id = id;
        this.title = title;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Book)) return false;
        Book book = (Book) o;
        return id == book.id && title.equals(book.title);
    }


    @Override
    public int hashCode() {
        return Objects.hash(id, title);
    }
//        return new HashCodeBuilder(17, 31). // two randomly chosen prime numbers
//                // if deriving: appendSuper(super.hashCode()).
//                        append(id).
//                append(title).
//                toHashCode();
//    }
    public int getId () {
        return id;
    }
    public String getTitle () {
        return title;
    }
}


