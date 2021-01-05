(define (caar x) (car (car x)))
(define (cadr x) (car (cdr x)))
(define (cdar x) (cdr (car x)))
(define (cddr x) (cdr (cdr x)))

; Some utility functions that you may find useful to implement.
(define (map proc items)
  (cond
    ((equal? items nil) nil)
    (else (cons (proc (car items)) (map proc (cdr items))))
    )
)

(define (helper s i)
    (cond
      ((equal? s nil) nil)
      (else (cons (list i (car s)) (helper (cdr s) (+ i 1))))
    )
)

(define (cons-all first rests)
  (define (helper elem_rest)
    (cons first elem_rest)
    )

  (map helper rests)
)

(define (zip pairs)
  (let ((first (map (lambda (pair) (car pair)) pairs))
        (rest (map (lambda (pair) (cadr pair)) pairs)))
    (list first rest))
)

;; Problem 17
;; Returns a list of two-element lists
(define (enumerate s)
  ; BEGIN PROBLEM 17
  (define i 0)
  (define (helper s i)
    (cond
      ((equal? s nil) nil)
      (else (cons (list i (car s)) (helper (cdr s) (+ i 1))))
      )
    )
  (helper s i)
  )
  ; END PROBLEM 17

;; Problem 18
;; List all ways to make change for TOTAL with DENOMS
(define (list-change total denoms)
  ; BEGIN PROBLEM 18
  (cond ((= total 0) (cons nil nil))
        ((or (< total 0) (null? denoms)) nil)
        ((< total (car denoms)) (list-change total (cdr denoms)))
        (else (append 
                (cons-all (car denoms) (list-change (- total (car denoms)) denoms))
                (list-change total (cdr denoms))
              )
        )
  )
)
  
  ; END PROBLEM 18

;; Problem 19
;; Returns a function that checks if an expression is the special form FORM
(define (check-special form)
  (lambda (expr) (equal? form (car expr))))

(define lambda? (check-special 'lambda))
(define define? (check-special 'define))
(define quoted? (check-special 'quote))
(define let?    (check-special 'let))

;; Converts all let special forms in EXPR into equivalent forms using lambda
(define (let-to-lambda expr)
  (cond ((atom? expr)
         ; BEGIN PROBLEM 19
         expr
         ; END PROBLEM 19
        )


        ((quoted? expr)
         ; BEGIN PROBLEM 19
         expr
         ; END PROBLEM 19
        )


        ((or (lambda? expr)
             (define? expr))
          (let ((form   (car expr))
               (params (cadr expr))
               (body   (cddr expr)))
           ; BEGIN PROBLEM 19
           (cons form (cons params (map let-to-lambda (cddr expr))))
           ; END PROBLEM 19
          )
        )



        ((let? expr)
          (let ((values (cadr expr))
               (body   (cddr expr)))
           ; BEGIN PROBLEM 19
            (let ((zipped_values (zip values)))
              (let  ((para_values (map let-to-lambda (cadr zipped_values)))
                    (func (append (list 'lambda (car zipped_values))
                            (map let-to-lambda body)
                          )
                    )
                    )
                    (cons func para_values)
              )
            )
           ; END PROBLEM 19
          )
        )



        (else
         ; BEGIN PROBLEM 19
          (map let-to-lambda expr)
         ; END PROBLEM 19
        )
  )
)
